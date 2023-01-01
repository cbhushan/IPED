% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function data_out = EPI_correct_slice_invLeastSq_image(dwi_slice, phase_dir, deltaB0_slice, echo_space, ...
                                             deform_scale_Factor, alpha, isfull, num_diffusion_dir, ...
                                             dwi_init, mask_undistorted, mask_dwi_slice)
%
% Estimates the diffusion-weighted images for optimal fit with only spatial
% regularization to the observed data. Operates one slice at a time, hence spatial
% regularization is only in-slice. 
% This works with at most four phase encode directions.
%
% dwi_slice - 3D matrix, say with size X*Y*n, where 3rd dimension (n) is the
%             different diffusion encoded images of the same slice 
%
% phase_dir - vector of integers (same length as 3rd dim of dwi_slice),
%                 representing the phase encode direction of nth diffusion
%                 encoded image in dwi_slice. Directions are represented
%                 as:
%                     1 - along first dimension of dwi_slice
%                     2 - along first dimension but opposite of 1 (replaces
%                         deltaB0 with -1*deltaB0)
%                     3 - along 2nd dimension of dwi_slice
%                     4 - along 2nd dimension but opposite of 3 
%
% deltaB0_slice - in Hz; 2D volume with size X*Y
%
% echo_space - echo spacing in sec
%
% deform_scale_Factor - defromation applied after upscaling the images by
%                       this factor. also check EPI_distort_fieldmap_image.m 
%
% alpha - scalar weights for spatial regularization
%
% isfull - Set to true if each EPI encoded images has been acquired with k>1 
%          number of phase encoding directions. If true, dwi_slices must have all 
%          EPI images with same phase encoding stacked together followed by 
%          images with other phase encoding. 3rd dim of dwi_slice should be 
%          n = k*d, where k is the number of phase encoding directions and d is the 
%          number of diffusion encodings (num_diffusion_dir). When n and d
%          are same (or k =1) then it is same as isfull set to false. 
% 
% num_diffusion_dir - Total number of diffusion encodings. Ignored when isfull is
%                     false. 
%
% mask_undistorted - (optional) Diffusion images are estimated only for
%                    pixels inside this mask (2D matrix, same size as
%                    deltaB0_slice) 
%
% mask_dwi_slice - (optional) data for pixels (from dwi_slice) outside this
%                   mask is not used for any estimation. (2D matrix, same
%                   size as deltaB0_slice)
%
%                                                      


% check if phase_dir is sane
if max(phase_dir)>4 || min(phase_dir)<1 || max(rem(phase_dir,1))>0
   error('phase_dir must be vector of integers in range [1 4]');
end


dwiSize = size(dwi_slice);
if length(dwiSize) == 3
   nVoxels = prod(dwiSize(1:2));
   ndim3 = dwiSize(end);
elseif length(dwiSize) == 2
   nVoxels = prod(dwiSize);
   ndim3 = 1;
else
   error('Number of dimension of dwi_slice must be either 2 or 3.')   
end


% check if isfull & related inputs are sane
if isfull && (ndim3==num_diffusion_dir)
   isfull = false; 
elseif isfull && mod(ndim3, num_diffusion_dir)~=0
   error('3rd dim of dwi_slice should be integer multiple of num_diffusion_dir')
elseif ~isfull
   num_diffusion_dir = ndim3;
end


% set masks
if exist('mask_dwi_slice', 'var')
   mask_dwi_slice = mask_dwi_slice>0;
else
   mask_dwi_slice = true(dwiSize(1:2));
end

if exist('mask_undistorted', 'var')
   mask_undistorted = mask_undistorted>0;
else
   mask_undistorted = true(dwiSize(1:2));
end

if ~exist('dwi_init', 'var')
   dwi_init = zeros([dwiSize(1:2) num_diffusion_dir]);
end
x_init = double(dwi_init(mask_undistorted(1:end, 1:end, ones(1,num_diffusion_dir))));
clear dwi_init

% check for all zero masks
if ~any(mask_undistorted(:)) || ~any(mask_dwi_slice(:)) 
   data_out = zeros([dwiSize(1:2) num_diffusion_dir]);
   return
end


% get different deformation matrices
if any(phase_dir==1)
   [~, deform_trsfrm1] = EPI_distort_fieldmap_image(dwi_slice(:,:,1), deltaB0_slice, echo_space, 1, ...
      deform_scale_Factor, mask_undistorted, mask_dwi_slice);
else
   deform_trsfrm1 = [];
end

if any(phase_dir==2)
   [~, deform_trsfrm2] = EPI_distort_fieldmap_image(dwi_slice(:,:,1), -1*deltaB0_slice, echo_space, 1, ...
      deform_scale_Factor, mask_undistorted, mask_dwi_slice);
else
   deform_trsfrm2 = [];
end

if any(phase_dir==3)
   [~, deform_trsfrm3] = EPI_distort_fieldmap_image(dwi_slice(:,:,1), deltaB0_slice, echo_space, 2, ...
      deform_scale_Factor, mask_undistorted, mask_dwi_slice);
else
   deform_trsfrm3 = [];
end

if any(phase_dir==4)
   [~, deform_trsfrm4] = EPI_distort_fieldmap_image(dwi_slice(:,:,1), -1*deltaB0_slice, echo_space, 2, ...
      deform_scale_Factor, mask_undistorted, mask_dwi_slice);
else
   deform_trsfrm4 = [];
end

dfSize = max([size(deform_trsfrm1); size(deform_trsfrm2); size(deform_trsfrm3); size(deform_trsfrm4)], [], 1);
if isempty(deform_trsfrm1)
   deform_trsfrm1 = sparse(dfSize(1), dfSize(2));
end
if isempty(deform_trsfrm2)
   deform_trsfrm2 = sparse(dfSize(1), dfSize(2));
end
if isempty(deform_trsfrm3)
   deform_trsfrm3 = sparse(dfSize(1), dfSize(2));
end
if isempty(deform_trsfrm4)
   deform_trsfrm4 = sparse(dfSize(1), dfSize(2));
end


% 1st order finite difference matrix 
finite_diff = createDWithPeriodicBoundary(dwiSize(1), dwiSize(2));
mask_ind = find(mask_undistorted<=0);
finite_diff(:, mask_ind) = [];
finite_diff([mask_ind; mask_ind+nVoxels],:) = [];


% create reshape index for two PE directions
temp = 1:length(phase_dir);
dir1_index = temp(phase_dir(:)==1);
dir2_index = temp(phase_dir(:)==2);
dir3_index = temp(phase_dir(:)==3);
dir4_index = temp(phase_dir(:)==4);


% create vector for observed data + zeros for regularization
d_len = ndim3*dfSize(1) + size(finite_diff, 1)*num_diffusion_dir;
d = zeros(d_len,1);
d(1:(ndim3*sum(mask_dwi_slice(:)))) = dwi_slice(mask_dwi_slice(1:end, 1:end, ones(1,ndim3)));


Afun = @(x, transp_flag)transform_X_invLeastSq_image(x, transp_flag, ...
                                    deform_trsfrm1, deform_trsfrm2, deform_trsfrm3, deform_trsfrm4,...
                                    alpha*finite_diff, isfull, num_diffusion_dir, ...
                                    dir1_index, dir2_index, dir3_index, dir4_index);                                 

[S, flag, relres] = lsqr(Afun, d, 1e-6, 1000, [], [], x_init);


% reshape estimated data back to original dimension
data_out = zeros([dwiSize(1:2) num_diffusion_dir]);
data_out(mask_undistorted(1:end, 1:end, ones(1,num_diffusion_dir))) = S;

end


function y = transform_X_invLeastSq_image(x, transp_flag, ...
                         deform_trsfrm1, deform_trsfrm2, deform_trsfrm3, deform_trsfrm4,...
                         finite_diff, isfull, num_diffusion_dir, ...
                         dir1_index, dir2_index, dir3_index, dir4_index)

if strcmp(transp_flag,'transp')      % y = A'*x
   
   if isfull
      FD_length = size(finite_diff, 1)*num_diffusion_dir;
   else
      FD_length = size(finite_diff, 1) * ...
                  (length(dir1_index) + length(dir2_index) + length(dir3_index) + length(dir4_index));
   end
   
   dwi_length = length(x) - FD_length;
      
   x_FD = reshape(x(dwi_length+1:end), 2*size(deform_trsfrm1, 2), []);
   y_FD = finite_diff' * x_FD;

   dwi = reshape(x(1:dwi_length), size(deform_trsfrm1, 1), []);
   y1 = deform_trsfrm1' * dwi(1:end, dir1_index);
   y2 = deform_trsfrm2' * dwi(1:end, dir2_index);
   y3 = deform_trsfrm3' * dwi(1:end, dir3_index);
   y4 = deform_trsfrm4' * dwi(1:end, dir4_index);
   
   if isfull
      if isempty(dir1_index), y1 = zeros(size(y1,1), num_diffusion_dir); end
      if isempty(dir2_index), y2 = zeros(size(y2,1), num_diffusion_dir); end
      if isempty(dir3_index), y3 = zeros(size(y3,1), num_diffusion_dir); end
      if isempty(dir4_index), y4 = zeros(size(y4,1), num_diffusion_dir); end
      
      y_deform_out = y1+y2+y3+y4;
   else
      y_deform = [y1 y2 y3 y4];
      clear y1 y2 y3 y4
      
      y_deform_out = zeros(size(y_deform));
      y_deform_out(1:end, [dir1_index dir2_index dir3_index dir4_index]) = y_deform;
      clear y_deform
   end
   y = y_deform_out + y_FD;
   y = y(:);
   
elseif strcmp(transp_flag,'notransp') % y = A*x
   dwi = reshape(x, size(deform_trsfrm1, 2), []);

   if isfull
      y1 = deform_trsfrm1 * dwi;
      y2 = deform_trsfrm2 * dwi;
      y3 = deform_trsfrm3 * dwi;
      y4 = deform_trsfrm4 * dwi;
      
      if isempty(dir1_index), y1 = zeros(size(y1,1), 0); end
      if isempty(dir2_index), y2 = zeros(size(y2,1), 0); end
      if isempty(dir3_index), y3 = zeros(size(y3,1), 0); end
      if isempty(dir4_index), y4 = zeros(size(y4,1), 0); end
   else
      y1 = deform_trsfrm1 * dwi(:, dir1_index);
      y2 = deform_trsfrm2 * dwi(:, dir2_index);
      y3 = deform_trsfrm3 * dwi(:, dir3_index);
      y4 = deform_trsfrm4 * dwi(:, dir4_index);
   end

   y_deform = [y1 y2 y3 y4];
   clear y1 y2 y3 y4
   y_deform_out = zeros(size(y_deform));
   y_deform_out(1:end, [dir1_index dir2_index dir3_index dir4_index]) = y_deform;
   clear y_deform
   
   % finite difference
   y_FD = finite_diff * dwi;
   
   y = [y_deform_out(:); y_FD(:)];
end

end
