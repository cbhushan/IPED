% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [mask, threshold, masked_data] = mask_head_pseudo(data_file, output_mask_file, if_aggresive)
% This function tries to detect head (including skull and skin) just based on
% intensity in a 3D volume. This function assumes that the head (or area of interest) and
% background have significant difference in their intensity values.  
%
% output_mask_file - (optional) file name of output mask
% if_aggresive - (optional) When true, it tries to be more aggressive in masking.
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



% Structuring elements
pix = 1;
[x1,y1,z1] = ndgrid(-pix:pix);
se1 = (sqrt(x1.^2 + y1.^2 + z1.^2) <=pix);

pix = 2;
[x1,y1,z1] = ndgrid(-pix:pix);
se2 = (sqrt(x1.^2 + y1.^2 + z1.^2) <=pix);

pix = 3;
[x1,y1,z1] = ndgrid(-pix:pix);
se3 = (sqrt(x1.^2 + y1.^2 + z1.^2) <=pix);
pix = 4;
[x1,y1,z1] = ndgrid(-pix:pix);
se4 = (sqrt(x1.^2 + y1.^2 + z1.^2) <=pix);

pix = 8;
[x1,y1,z1] = ndgrid(-pix:pix);
se8 = (sqrt(x1.^2 + y1.^2 + z1.^2) <=pix);

data_orig = data;
data_size = size(data.img);
[data.img, lw, hg] = normalize_intensity(double(data.img), [5 95]);

% smooth the image to get robust threshold
opt = struct('kappa', 15, 'lambda', 1/8, 'eq', 1);
data_smooth = anisotropic_diffusion_filter3d(double(data.img), 4, opt)*1200;

% find a reasonable threshold based on first derivative of sorted values
I_s = sort(data_smooth(:), 'ascend');
d = floor(length(I_s)/5000);  % pick 5000 equidistant samples
x = smooth_moving(double(I_s(1:d:end)), 31);
x_diff = diff(x);
threshold_smooth = find_thrshold_MASK_HEAD_PSEUDO(x, x_diff, 0.85);
mask_smooth = data_smooth>threshold_smooth;

% padd zeros to avoid weird behaviour at boundaries
pad_size = 4;
mask_smooth_pad = padarray(mask_smooth, [1 1 1]*pad_size);

% remove isolated pixels
msk_tmp = imerode(mask_smooth_pad>0, se3)>0;
msk_tmp = largest_connected_component_MASK_HEAD_PSEUDO(msk_tmp, 6);
msk_tmp = imdilate(msk_tmp>0, se3)>0;
msk_tmp = msk_tmp & mask_smooth_pad>0;
msk_tmp = imfill(msk_tmp>0, 'holes');

% Try to get back small thin pieces lost during imerode
msk_tmp = imdilate(msk_tmp>0, se2)>0 & (mask_smooth_pad>0);
msk_tmp = largest_connected_component_MASK_HEAD_PSEUDO(msk_tmp, 6);
mask_smooth = msk_tmp(pad_size+1:data_size(1)+pad_size, pad_size+1:data_size(2)+pad_size, pad_size+1:data_size(3)+pad_size)>0;
clear msk_tmp mask_smooth_pad

mask = data;
if exist('if_aggresive', 'var') && if_aggresive
   I_s = sort(data.img(mask_smooth), 'ascend');
   threshold = max(I_s(floor(length(I_s)*0.055)), threshold_smooth/1200);
   mask.img = data.img>=threshold;
   
   % padd zeros to avoid weird behaviour at boundaries
   mask_pad = padarray(mask.img, [1 1 1]*pad_size);
   mask_pad_fill = imfill(mask_pad>0, 'holes');
   
   % remove isolated pixels
   msk_tmp = imerode(mask_pad>0, se3)>0;
   msk_tmp = largest_connected_component_MASK_HEAD_PSEUDO(msk_tmp, 6);
   msk_tmp = imdilate(msk_tmp>0, se4)>0;
   msk_tmp = msk_tmp & mask_pad>0;
   msk_tmp = imfill(msk_tmp>0, 'holes');
   
   % check if we are missing too many voxels - to avoid cases when aggressive thresholding only
   % gets ventricles in T2-weighted images
   [bb_1, bb_2, bb_3] = find_bounding_box(msk_tmp);
   vol_ratio = sum(msk_tmp(:))/((bb_1(end)-bb_1(1))*(bb_2(end)-bb_2(1))*(bb_3(end)-bb_3(1)));
   if vol_ratio<0.3
      threshold = find_thrshold_MASK_HEAD_PSEUDO(x, x_diff, 0.65); % relax a little bit
      mask.img = data_orig.img>=(threshold/1200);
      
      % padd zeros to avoid weird behaviour at boundaries
      mask_pad = padarray(mask.img, [1 1 1]*pad_size);
      
      % remove isolated pixels
      msk_tmp = imerode(mask_pad>0, se2)>0;
      msk_tmp = largest_connected_component_MASK_HEAD_PSEUDO(msk_tmp, 6);
      msk_tmp = imdilate(msk_tmp>0, se2)>0;
      msk_tmp = msk_tmp & mask_pad>0;
   end
   
   
   % close holes only deep inside
   msk_tmp = imdilate(msk_tmp>0, se1) & mask_pad_fill>0;
   msk_tmp = imfill(msk_tmp>0, 'holes');
   msk_tmp = imdilate(msk_tmp>0, se1) & mask_pad_fill>0;
   msk_tmp = imfill(msk_tmp>0, 'holes');
   
   mask.img =  msk_tmp(pad_size+1:data_size(1)+pad_size, pad_size+1:data_size(2)+pad_size, pad_size+1:data_size(3)+pad_size);
   clear msk_tmp mask_pad_fill mask_pad
else
   mask.img = mask_smooth;
   threshold = threshold_smooth/1200;
end

% make BrainSuite compatible
mask.img = uint8(mask.img>0)*255;

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
function out = largest_connected_component_MASK_HEAD_PSEUDO(mask, conn)

cc = bwconncomp(mask, conn);
cc_size = [];
for k = 1:cc.NumObjects
   cc_size(k) = length(cc.PixelIdxList{k});
end
[~,IX] = sort(cc_size,'descend');
out = zeros(size(mask));
out(cc.PixelIdxList{IX(1)}) = 1;

end


function threshold = find_thrshold_MASK_HEAD_PSEUDO(x, x_diff, const)

ind = find(x_diff(500:end)>const, 1, 'first'); % ignore lower 10% intensities (wierd -ve/0 behaviour)

if ~isempty(ind)
   threshold = x(ind+499);
   
else % flat 10th percentile as threshold
   threshold = x(floor(length(x)/10));
end
end

function smooth_data = smooth_moving(data, windowSize)

if mod(windowSize,2)==0 % even windowSize
   windowSize = windowSize+1;
end

if ndims(data)~=2
   error('data must be vector.')
end

smooth_data = filter(ones(1,windowSize)/windowSize, 1, data);

end
