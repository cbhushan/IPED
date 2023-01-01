% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [data_out, deform_transform, i, j, s, m, n] = EPI_distort_fieldmap_image(data, deltaB0, echo_space, distort_dim,...
                                                                                   deform_scale_factor, data_mask, mask_distorted)
% Distorts the data using field-map in image domain (by expressing the
% distortion as linear system). 
%  
%  Could be VERY memory intensive.
% 
%   data - 2D volume, with 1st dimension as phase encode &
%            second dimension as readout direction (One slice at a time)
%   deltaB0 - in Hz; 2D volume (dimensions same as DATA)
%   echo_spacing  - in sec
%   distort_dim - 1 or 2; Dimension along which the input data should be
%                 distorted (for opposite distortion in the same dimension
%                 pass -1*deltaB0 as input)                 
%   deform_scale_factor - Initial distortion map is computed by upscaling by
%                 this factor - helps in smoother & realistic distortion by 
%                 using a linear approximation for modeling inhomogeneity in 
%                 each voxel. Use power of 2 for better results.
%   data_mask - (optional) area outside the mask will be discarded
%   mask_distorted - (optional) area outside this mask will be always set
%                    to zero in data_out. 
%
% Implementation notes: 
%    deform_transform upsamples the input image (using bilinear interpolation), and then
%    computes mapping in between upsampled images. However, the mapping is then downsampled to
%    original voxel resolution by averaging the effect over each voxel (at original resolution).
%    This helps in mixing of signal over multiple voxels. This is similar to using linear 
%    approximation for modeling inhomogeneity in each voxel. 
%    
%    When using deform_transform repeatedly, like in iterative alogrithms, it may be useful to use
%    masks (data_mask and mask_distorted) which can reduce memory requirement and run time
%    substantially. 
%

sz = size(data);
num_vox = sz(1)*sz(2);


if distort_dim<1 || distort_dim>2
   error('distort_dim must be either 1 or 2')
   
elseif distort_dim==1
   shift = double(deltaB0) * sz(1) * double(echo_space);
   
   [x_g, y_g] = ndgrid(1:sz(1), 1:sz(2));
   x_new = x_g + shift;
   
   % upscale the map to get more contineous map
   [interp_mat, sz_up] = interpolation_matrix_2D(sz, [deform_scale_factor 1]);
   
   x_g_up = interp_mat * double(x_g(:));
   x_g_up = reshape(x_g_up, sz_up);
   
   x_new_up = interp_mat * double(x_new(:));
   x_new_up = reshape(x_new_up, sz_up);
   
   y_g_up = interp_mat * double(y_g(:));
   y_g_up = reshape(y_g_up, sz_up);

   % weights for linear interpolation
   new_floor_weight = ceil(x_new_up)-x_new_up;
   new_ceil_weight = 1-new_floor_weight;

   % Find voxel to voxel (index to index) map and corresponding weights
   x_new_ceil = ceil(x_new_up);
   ceil_invalid_mask = x_new_ceil>sz(1) | x_new_ceil<1;
   ceil_new_vox_ind = sub2ind(sz, x_new_ceil(~ceil_invalid_mask),  y_g_up(~ceil_invalid_mask));
   ceil_vox_ind = sub2ind(sz_up, x_g_up(~ceil_invalid_mask).*deform_scale_factor-(deform_scale_factor-1),  y_g_up(~ceil_invalid_mask));
   new_ceil_weight = new_ceil_weight(~ceil_invalid_mask); % remove points outside the grid
   
   x_new_floor = floor(x_new_up);
   floor_invalid_mask = x_new_floor>sz(1) | x_new_floor<1;
   floor_new_vox_ind = sub2ind(sz, x_new_floor(~floor_invalid_mask),  y_g_up(~floor_invalid_mask));
   floor_vox_ind = sub2ind(sz_up, x_g_up(~floor_invalid_mask).*deform_scale_factor-(deform_scale_factor-1),  y_g_up(~floor_invalid_mask));
   new_floor_weight = new_floor_weight(~floor_invalid_mask); % remove points outside the grid
   
elseif distort_dim==2
   shift = double(deltaB0) * sz(2) * double(echo_space);
   
   [x_g, y_g] = ndgrid(1:sz(1), 1:sz(2));
   y_new = y_g + shift;
   
   % upscale the map to get more contineous map
   [interp_mat, sz_up] = interpolation_matrix_2D(sz, [1 deform_scale_factor]);
   
   x_g_up = interp_mat * double(x_g(:));
   x_g_up = reshape(x_g_up, sz_up);
   
   y_new_up = interp_mat * double(y_new(:));
   y_new_up = reshape(y_new_up, sz_up);
   
   y_g_up = interp_mat * double(y_g(:));
   y_g_up = reshape(y_g_up, sz_up);

   % weights for linear interpolation
   new_floor_weight = ceil(y_new_up)-y_new_up;
   new_ceil_weight = 1-new_floor_weight;

   % Find voxel to voxel (index to index) map and corresponding weights
   y_new_ceil = ceil(y_new_up);
   ceil_invalid_mask = y_new_ceil>sz(2) | y_new_ceil<1;
   ceil_new_vox_ind = sub2ind(sz, x_g_up(~ceil_invalid_mask),  y_new_ceil(~ceil_invalid_mask));
   ceil_vox_ind = sub2ind(sz_up, x_g_up(~ceil_invalid_mask),  y_g_up(~ceil_invalid_mask).*deform_scale_factor-(deform_scale_factor-1));
   new_ceil_weight = new_ceil_weight(~ceil_invalid_mask); % remove points outside the grid
   
   y_new_floor = floor(y_new_up);
   floor_invalid_mask = y_new_floor>sz(2) | y_new_floor<1;
   floor_new_vox_ind = sub2ind(sz, x_g_up(~floor_invalid_mask),  y_new_floor(~floor_invalid_mask));
   floor_vox_ind = sub2ind(sz_up, x_g_up(~floor_invalid_mask),  y_g_up(~floor_invalid_mask).*deform_scale_factor-(deform_scale_factor-1));
   new_floor_weight = new_floor_weight(~floor_invalid_mask); % remove points outside the grid
end

% create sparse matrix 
i = [ceil_new_vox_ind; floor_new_vox_ind];
j = [ceil_vox_ind; floor_vox_ind];
s = [new_ceil_weight; new_floor_weight]/(deform_scale_factor);
m = num_vox;
n = prod(sz_up);

deform_transform = sparse(i,j,s,m,n);

% include upsampling matrix in the deformation matrix 
deform_transform = deform_transform * interp_mat; 

% throw away some columns
if exist('data_mask', 'var')
   mask_ind = find(data_mask<=0);
   deform_transform(:, mask_ind) = [];
   data(mask_ind) = [];
end

% throw away some rows
if exist('mask_distorted', 'var')
   mask_ind = find(mask_distorted<=0);
   deform_transform(mask_ind,:) = [];
end


% distort the data & then reshape output properly
d_trans = deform_transform*double(data(:));

data_out = zeros(sz);
if exist('mask_distorted', 'var')
   data_out_mask = mask_distorted>0;
else
   data_out_mask = true(sz);
end

data_out(data_out_mask) = d_trans;

end
