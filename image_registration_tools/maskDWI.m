% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [mask_file, mask_less_csf_file, mask_bfc] = maskDWI(DWI, outfile_base, varargin)
% Tries to estimate mask from 4D DWI volume. DWI should NOT be (severely) affected by motion. 
% B0 drift and severe bias-field/non-uniformity can also affect quality of result. 
%
% Usage: 
%   [mask, mask_less_csf, mask_bfc] = maskDWI(DWI, out_file_base, bmat_file)
%   [mask, mask_less_csf, mask_bfc] = maskDWI(DWI, out_file_base, bmat_file, masking_approach)
%   [mask, mask_less_csf, mask_bfc] = maskDWI(DWI, out_file_base, bvec_file, bval_file)
%   [mask, mask_less_csf, mask_bfc] = maskDWI(DWI, out_file_base, bvec_file, bval_file, masking_approach)
%
% masking_approach is a numeric flag which defines the masking approach: 
%   1  -  'intensity'; Heuristic approach based on intensities of b=0 image with mask_head_pseudo(). This is
%         the oldest method implemented in BDP - still works best for some datasets.
%   2  -  'hist'; Using peaks and valleys of the intesity histogram of (b0 + meanDWI) image. This makes
%         direct intuitive sense, but fails when data is already partially masked, like that in
%         HCP data.
%   3  -  (Not implemented) Spectral clustering based segmentation should be robust to a lot of
%         issues?


nbins = 300; 
parzen_width = 40;
nthreads = 6;
valley_slp_thresh = 20; 
masking_approach = 2; % use histogram based masking as default

% check input options
if nargin==3
   de_out = checkDiffusionEncodingScheme(varargin{1});
elseif nargin==4
   if isscalar(varargin{2}) && ~ischar(varargin{2})
      masking_approach = varargin{2};
      de_out = checkDiffusionEncodingScheme(varargin{1});
   else
      de_out = checkDiffusionEncodingScheme(varargin{1}, varargin{2});
   end
elseif nargin==5
   de_out = checkDiffusionEncodingScheme(varargin{1}, varargin{2});
   masking_approach = varargin{3};
else
   error('Incorrect number of inputs. Check usage.')
end

% load input image
if ischar(DWI)
   data = load_untouch_nii_gz(DWI, true);
else
   if ~isfield(DWI,'untouch') || DWI.untouch ~= 1
      error('Please use ''load_untouch_nii.m'' to load the file.');
   end
   data = DWI;
   clear DWI
end

outfile_base = remove_extension(outfile_base);

% separate b=0 and DWIs
if ndims(data.img)==4
   b0 = data;
   b0.hdr.dime.dim(1) = 3;
   b0.hdr.dime.dim(5) = 1;
   dwi_mean = b0;
   dwi_mean.img = double(mean(data.img(:,:,:,~de_out.zero_bval_mask), 4));
   
   % stick to 1st b=0 instead of mean b=0 as 1st one will be used for registration 
   b0.img = double(data.img(:,:,:,find(de_out.zero_bval_mask,1)));
   
elseif ndims(data.img)==3
   fprintf('\n3D DWI volume found! BDP will assume that it is b=0 image.\n');
   b0 = data;
   dwi_mean = [];
else
   error('DWI is not 3D or 4D volume - not supported.')
end
clear data

if ~isempty(b0.img) && ~isempty(dwi_mean)
   switch masking_approach
      case 1
         mask = mask_head_pseudo(b0, true);
         
      case 2
         b0.img = normalize_intensity(b0.img, [3 97]);
         dwi_mean.img = normalize_intensity(dwi_mean.img, [3 97]);
         temp = b0;
         temp.img = (b0.img + dwi_mean.img)/2;
         mask = maskHeadPseudoHist(temp, true);
      otherwise
         error('Unknown/unimplemented masking_approach: %d', masking_approach);
   end
   [mask_dwi, dwi_thresh] = maskHeadPseudoHist(dwi_mean, true); % aggresive DWI masking
   
   mask_less_csf = mask_dwi;
   mask_less_csf.img = (imerode(mask.img>0, strel_sphere(3)) | mask_dwi.img>0) & mask.img>0;       
   mask_less_csf.img = largest_connected_component(mask_less_csf.img, 6);
   
   if nargout>2  % try to get bias-corrected estimate     
      sure_mask = imerode(mask.img>0 & mask_dwi.img>0, strel_sphere(1));
      dwi_bfc.img = dwi_mean.img ./ b0.img;    % always <=1, unless noise
      dwi_bfc.img(~isfinite(dwi_bfc.img)) = 1; % set to 1 to enhance contrast b/w brain and background
      
      % background could be super noisy - get rid of background
      
      % make brain part smoother
      opt = struct('kappa', 12, 'lambda', 1/8, 'eq', 1, 'voxRes', b0.hdr.dime.pixdim(2:4));
      dwi_bfc_smooth = anisotropic_diffusion_filter3d(double(dwi_bfc.img), 2, opt);
      dataN = normalize_intensity(dwi_bfc_smooth, [0 99.9]);
      
      % find valley in histogram
      [~,~, hist_data, bin_cntr] = histcParzen(dataN(:), [0 1], nbins, parzen_width, nthreads);
      gh = gradient(hist_data); % +ve 2nd derivative = Valley
      [bin0cross, bin0cross_slp] = zeroCrossing(gh, bin_cntr);
      ind = find(bin0cross_slp>valley_slp_thresh);
      
      if isempty(ind)
         mask_bfc_FG = sure_mask;
      else
         threshold = bin0cross(ind(1)); % 1st large zero crossing with +ve slope
         mask_bfc_FG = dataN<threshold;
         mask_bfc_FG = mask_bfc_FG | sure_mask;
         mask_bfc_FG = imfill(mask_bfc_FG, 'holes'); % close holes
         
         % cleanup mask - remove isolated components
         data_size = size(mask_bfc_FG);
         
         % padd zeros to avoid weird behaviour at boundaries
         pad_size = 4;
         mask_data_pad = padarray(mask_bfc_FG, [1 1 1]*pad_size)>0;
         
         % remove isolated pixels
         msk_tmp = imerode(mask_data_pad, strel_sphere(2));
         msk_tmp = largest_connected_component(msk_tmp, 6);
         msk_tmp = imdilate(msk_tmp, strel_sphere(3));
         msk_tmp = msk_tmp & mask_data_pad;
         msk_tmp = largest_connected_component(msk_tmp, 6);
         msk_tmp = imfill(msk_tmp>0, 'holes');
         
         mask_bfc_FG = msk_tmp(pad_size+1:data_size(1)+pad_size, pad_size+1:data_size(2)+pad_size, pad_size+1:data_size(3)+pad_size)>0;
         clear msk_tmp mask_data_pad
      end
      
      % mask_FG will have some noisy voxels - heuristic removal
      opt = struct('kappa', 12, 'lambda', 1/8, 'eq', 1, 'voxRes', b0.hdr.dime.pixdim(2:4));
      dwi_bfc_smooth = anisotropic_diffusion_filter3d(double(dwi_bfc.img), 5, opt);
      
      [hy, hx, hz] = gradient(dwi_bfc_smooth, b0.hdr.dime.pixdim(2), b0.hdr.dime.pixdim(2), b0.hdr.dime.pixdim(2));
      tt_g = sqrt(hy.^2 + hx.^2 + hz.^2);
      g_thresh = prcentileIntensity(tt_g(mask_bfc_FG), 90); % 90th percentile of gradient inside FG
      
      pad_size = 4;
      mask_data_pad = padarray(mask_bfc_FG, [1 1 1]*pad_size)>0;
      temp_g = padarray(tt_g, [1 1 1]*pad_size);
      
      msk_tmp = imerode(mask_data_pad, strel_sphere(3));
      msk_tmp = largest_connected_component(msk_tmp, 6);
      msk_tmp = imdilate(msk_tmp, strel_sphere(1)) & (temp_g<g_thresh) & mask_data_pad;
      msk_tmp = imdilate(msk_tmp, strel_sphere(1)) & (temp_g<g_thresh) & mask_data_pad;
      msk_tmp = imdilate(msk_tmp, strel_sphere(1)) & (temp_g<g_thresh) & mask_data_pad;
      msk_tmp = largest_connected_component(msk_tmp, 6);
      
      mask_bfc = mask;
      mask_bfc.img = msk_tmp(pad_size+1:data_size(1)+pad_size, pad_size+1:data_size(2)+pad_size, pad_size+1:data_size(3)+pad_size)>0;
      mask_bfc.img = mask_bfc.img | sure_mask;
   end
   
elseif ~isempty(dwi_mean.img)
   [mask, dwi_thresh] = maskHeadPseudoHist(dwi_mean, true);
   mask_less_csf = mask;
   mask_bfc = mask;

elseif ~isempty(b0.img)
   [mask, b0_thresh] = maskHeadPseudoHist(b0, true);
   mask_less_csf = mask;
   mask_bfc = mask;
end

% write outputs
mask.hdr.dime.scl_slope = 0;
mask.img = uint8(mask.img>0)*255; % make BrainSuite compatible
mask_file = save_untouch_nii_gz(mask, [outfile_base '.mask.nii.gz'], 2);

mask_less_csf.hdr.dime.scl_slope = 0;
mask_less_csf.img = uint8(mask_less_csf.img>0)*255; % make BrainSuite compatible
mask_less_csf_file = save_untouch_nii_gz(mask_less_csf, [outfile_base '.less_csf.mask.nii.gz'], 2);

end


% Find the biggest connected component
function out = largest_connected_component(mask, conn)

cc = bwconncomp(mask, conn);
if cc.NumObjects>1
   cc_size = [];
   for k = 1:cc.NumObjects
      cc_size(k) = length(cc.PixelIdxList{k});
   end
   [~,IX] = sort(cc_size,'descend');
   out = false(size(mask));
   out(cc.PixelIdxList{IX(1)}) = true;
else
   out = mask;
end
end


function [bin0cross, bin0cross_slp, bin0kiss] = zeroCrossing(x, bin)

n = length(x);
ind0_kiss = [];

% First handle exact zeros
ind0 = find(x==0);
if ~isempty(ind0)
   tempind = ind0;
   if tempind(1)==1
      tempind(1) = [];
      ind0_kiss = [ind0_kiss 1];
   end
   if tempind(end)==n
      tempind(end) = [];
      ind0_kiss = [ind0_kiss n];
   end
   
   % following misses contiguous zero values of x, but thats OK here
   px = x(tempind+1) .* x(tempind-1);
   ind0_kiss = [ind0_kiss tempind(px>0)];
end

% look for zero crossings
ap = x(1:end-1) .* x(2:end); % adjacent points should have different signs
ind1 = find(ap<0);

frc = abs(x(ind1) ./ (x(ind1+1) - x(ind1)));
bin0cross = bin(ind1) + (frc .* (bin(ind1+1) - bin(ind1))); % actual location of zero crossing

% add exact zeros & sort
bin0cross = sort([bin0cross(:); vect(bin(ind0))], 'ascend');
bin0kiss = sort(vect(bin(ind0_kiss)), 'ascend');

% compute gradients 
dx = gradient(x, bin(2)-bin(1));
bin0cross_slp = interp1(bin, dx, bin0cross);

end

function t = prcentileIntensity(d, prct)
d = sort(d(:), 'ascend');
ind = ceil(prct/100*length(d));
t = d(ind);
end

function img = stretch_intensity(img, msk)
if nargin==1
   msk = true(size(img));
end

l = min(img(msk));
h = max(img(msk));
img = (img-l)/(h-l);
end

function h = log_kernel2D(p2, p3)
% Laplacian of Gaussian
% h = log_kernel(size, sigma)

% first calculate Gaussian
siz   = (p2-1)/2;
std2   = p3^2;

[x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
arg   = -(x.*x + y.*y)/(2*std2);

h     = exp(arg);
h(h<eps*max(h(:))) = 0;

sumh = sum(h(:));
if sumh ~= 0,
   h  = h/sumh;
end;
% now calculate Laplacian
h1 = h.*(x.*x + y.*y - 2*std2)/(std2^2);
h     = h1 - sum(h1(:))/prod(p2); % make the filter sum to zero

end

function h = log_kernel3D(siz)
% Laplacian of Gaussian
sig = siz/2/2.354;
siz   = (siz-1)/2;
[x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));
h = h/sum(h(:));
arg = (x.*x/sig(1)^4 + y.*y/sig(2)^4 + z.*z/sig(3)^4 - ...
   (1/sig(1)^2 + 1/sig(2)^2 + 1/sig(3)^2));
h = arg.*h;
h = h-sum(h(:))/prod(2*siz+1);
end
