% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [mask, threshold, masked_data] = maskHeadPseudoOtsu(data_file, output_mask_file, if_aggresive)
% This function tries to detect head (including skull and skin) just based
% on intensity in a 3D volume. This function uses Otsu's method and
% implicitly assumes that the head (or area of interest) and background
% have significant differences in their intensity values.
%
% This function is drop-in replacement for mask_head_pseudo.m and should
% outperform it in most cases. Also see mask_head_pseudo.m and maskHeadPseudoHist.m
%
% Inputs: 
%    data_file - nifti struct or filename
%    output_mask_file - (optional) file name of output mask
%    if_aggresive - (optional) When true, produces more eroded mask.
%
% Outputs:
%    mask - nifti struct of mask
%    threshold - Approximate intensity threshold. See note below for
%                explantion.
%    masked_data - nifti struct of masked data
%
% Always check output - may NOT make any sense at all!
%

% check input options
if nargin==2 && ~ischar(output_mask_file)
   if_aggresive = output_mask_file;
   clear output_mask_file
end

% load input image and some sanity test
if ischar(data_file)
   data = load_untouch_nii_gz(data_file, true);
else
   if ~isfield(data_file,'untouch') || data_file.untouch ~= 1
      error('Please use ''load_untouch_nii.m'' to load the file.');
   end
   data = data_file;
   clear data_file
end

if ndims(data.img)>3
   error('data_file must be 3D volume.')
end

% check if the input is a mask itself
unq_intensities = unique(data.img(:));
if length(unq_intensities)<=2
   mask = data;
   
   if length(unq_intensities)==2
      threshold = min(unq_intensities);
   elseif length(unq_intensities)==1
      threshold = unq_intensities-1;
   end
   mask.img = mask.img>threshold;
   
   % make BrainSuite compatible
   mask.img = uint8(mask.img>0)*255;
   
   if exist('output_mask_file', 'var') 
      mask.hdr.dime.scl_slope = 0;
      save_untouch_nii_gz(mask, output_mask_file, 2);
   end
   
   if nargout>2
      masked_data = data;
      masked_data.img(mask.img<=0) = 0;
   end
   
   return;
end

% Structuring elements; highly anisotropic voxels may cause issues
se1 = strel_sphere(1);
se2 = strel_sphere(2);
se3 = strel_sphere(3);
se4 = strel_sphere(4);

data_orig = data;
data_size = size(data.img);
[data.img, lw, hg] = normalize_intensity(double(data.img), [5 95]);

% anisotropic smooth the image to get robust threshold - this helps to
% "tame" the histogram and results into better otsu threshold.
% Also helps a little bit with local image non-uniformity
voxResRelative = double(data.hdr.dime.pixdim(2:4))/min(double(data.hdr.dime.pixdim(2:4)));
opt = struct('kappa', 17, 'lambda', 1/8, 'eq', 1, 'voxRes', voxResRelative);
data_smooth = anisotropic_diffusion_filter3d(double(data.img), 4, opt);
data_smooth = normalize_intensity(data_smooth, [0 100]); % put data in range [0 1]

% OTSU threshold and its effectiveness
[threshold_level, em] = graythresh(data_smooth(:)); 

% em (effectiveness metric) is ~0.75 for uniformally random variable
% em is typically ~0.9 for nicely separated gaussian distribution
if em < 0.75
   hdr = 'Check (pseudo) masking';
   msg = ['It is highly likely that the pseudo masking performed is not reasonable. '...
      'Effectiveness metric of thresholding was ' num2str(em) ...
      '. Please check output and, if required, use manual masking!'];
   bdpPrintWarning(hdr, msg);
end
mask_smooth = data_smooth>threshold_level;

% Compute an approximate threshold
% As smoothed version of image is used for thresholding this would always
% be approximate! Max or min wont work because of smoothing!
tmp = sort(data_orig.img(mask_smooth<=0), 'ascend');
threshold = tmp(floor(0.99*length(tmp))); % 99th percentile of background

% padd zeros to avoid weird behaviour at boundaries
pad_size = 4;
mask_smooth_pad = padarray(mask_smooth, [1 1 1]*pad_size);

% remove isolated pixels
msk_tmp = imerode(mask_smooth_pad>0, se3)>0;
msk_tmp = largest_connected_component_maskHeadPseudoOtsu(msk_tmp, 6);
if exist('if_aggresive', 'var') && if_aggresive
   msk_tmp = imdilate(msk_tmp>0, se2)>0;
else
   msk_tmp = imdilate(msk_tmp>0, se4)>0;
end
msk_tmp = msk_tmp & mask_smooth_pad>0;
msk_tmp = imfill(msk_tmp>0, 'holes');
mask_smooth = msk_tmp(pad_size+1:data_size(1)+pad_size, pad_size+1:data_size(2)+pad_size, pad_size+1:data_size(3)+pad_size)>0;
clear msk_tmp mask_smooth_pad


mask = data;
mask.img = mask_smooth;
mask.img = uint8(mask.img>0)*255; % make BrainSuite compatible

if exist('output_mask_file', 'var')
   mask.hdr.dime.scl_slope = 0;
   save_untouch_nii_gz(mask, output_mask_file, 2);
end

if nargout>2
   masked_data = data_orig;
   masked_data.img(mask.img<=0) = 0;
end


end

% Find the biggest connected component
function out = largest_connected_component_maskHeadPseudoOtsu(mask, conn)

cc = bwconncomp(mask, conn);
cc_size = [];
for k = 1:cc.NumObjects
   cc_size(k) = length(cc.PixelIdxList{k});
end
[~,IX] = sort(cc_size,'descend');
out = zeros(size(mask));
out(cc.PixelIdxList{IX(1)}) = 1;

end
