% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [data_out, c_out] = EPI_correct_slice_invLeastSq_spharm(dwi_slice, phase_dir, deltaB0_slice, echo_space, ...
                                                       deform_scale_factor, bvec, shOrder, alpha, isfull, ...
                                                       num_diffusion_dir, options)
%
% Estimates the spharmonic coefficients (for diffusion images) with
% regularization & optimal fit to the observed data. Operates one slice at
% a time. Works with at most four orthogonal phase encode directions.
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
%                       this factor. also chec EPI_distort_fieldmap_image.m 
%
% alpha - vector of length two, which are weights for regularization. First 
%         element is regularization weight for imposing spherical smoothness. 
%         second weight is for imposing spatial smoothness.
%
% isfull - if true, the dwi_slices should have all images with same phase
%          encode stacked together followed by images with other phase
%          encodes. 3rd dim of dwi_slice should be n = k*d, where k is the
%          number of fully sampled phase encoding directions and d is the 
%          number of diffusion encodings (num_diffusion_dir). When n and d
%          are same (or k =1) then it is same as isfull set to false. 
% 
% num_diffusion_dir - number of diffusion encodings. Ignored when isfull is
%                     false. 
%
% options.mask_undistorted - (optional) Diffusion images are estimated only for
%                    pixels inside this mask (2D matrix, same size as
%                    deltaB0_slice) 
%
% options.mask_dwi_slice - (optional) data for pixels (from dwi_slice) outside this
%                   mask is not used for any estimation. (2D matrix, same
%                   size as deltaB0_slice)
%
% options.spherical_weight - 
%
%


% check if phase_dir is sane
if max(phase_dir)>4 || min(phase_dir)<1
   error('phase_dir must be integer in range [1 4]');
end

if max(rem(phase_dir,1))>0
   error('phase_dir must be an integer in range [1 4]');
end

dwiSize = size(dwi_slice);
nVoxels = prod(dwiSize(1:end-1));
ndim3 = dwiSize(end);


% check if isfull & related inputs are sane
if isfull && (ndim3==num_diffusion_dir)
   isfull = false; 
elseif isfull && mod(ndim3, num_diffusion_dir)~=0
   error('3rd dim of dwi_slice should be multiple of num_diffusion_dir')
elseif ~isfull
   num_diffusion_dir = ndim3;
end

dwi_slice = double(dwi_slice);
deltaB0_slice = double(deltaB0_slice);

% set masks
if exist('options', 'var') && isfield(options, 'mask_dwi_slice') && ~isempty(options.mask_dwi_slice)
   mask_dwi_slice = options.mask_dwi_slice(1:end, 1:end, ones(1,ndim3))>0;
else
   mask_dwi_slice = true(dwiSize);
end

if exist('options', 'var') && isfield(options, 'mask_undistorted') && ~isempty(options.mask_undistorted)
   mask_undistorted = options.mask_undistorted(1:end, 1:end, ones(1,ndim3))>0;
else
   mask_undistorted = true(dwiSize);
end

if exist('options', 'var') && isfield(options, 'spherical_weight') && ~isempty(options.spherical_weight)
   spherical_weight = double(options.spherical_weight);
else
   spherical_weight = ones(dwiSize(1:end-1));
end


% modified SH basis
[theta_DE, phi_DE, ~] = cart2sph(bvec(:,1),bvec(:,2),bvec(:,3));
Y_sh = spharm_Hardi(shOrder, theta_DE, phi_DE); % basis matrix
% if isfull
%    Y_sh = Y_sh(1:ndim3/2,:);
% end
Ysize = size(Y_sh);
clear theta_DE phi_DE


% check for all zero masks
if ~any(mask_undistorted(:)) || ~any(mask_dwi_slice(:)) 
   data_out = zeros([dwiSize(1:2) num_diffusion_dir]);
   c_out = zeros([dwiSize(1:2) Ysize(2)]);
   return
end


% get different deformation matrices
[~, deform_trsfrm1] = EPI_distort_fieldmap_image(dwi_slice(:,:,1), deltaB0_slice, echo_space, 1, ...
                                                  deform_scale_factor, mask_undistorted(:,:,1), mask_dwi_slice(:,:,1));

[~, deform_trsfrm2] = EPI_distort_fieldmap_image(dwi_slice(:,:,1), -1*deltaB0_slice, echo_space, 1, ...
                                                  deform_scale_factor, mask_undistorted(:,:,1), mask_dwi_slice(:,:,1));
                                               
[~, deform_trsfrm3] = EPI_distort_fieldmap_image(dwi_slice(:,:,1), deltaB0_slice, echo_space, 2, ...
                                                  deform_scale_factor, mask_undistorted(:,:,1), mask_dwi_slice(:,:,1));

[~, deform_trsfrm4] = EPI_distort_fieldmap_image(dwi_slice(:,:,1), -1*deltaB0_slice, echo_space, 2, ...
                                                  deform_scale_factor, mask_undistorted(:,:,1), mask_dwi_slice(:,:,1));
                                               
dfSize = size(deform_trsfrm1);

% [~, deform_trsfrm1] = EPI_distort_fieldmap_image(dwi_slice(:,:,1), deltaB0_slice, echo_space, ...
%                                                  deform_scale_factor, mask_undistorted(:,:,1), mask_dwi_slice(:,:,1));
%                                               
% [~, deform_trsfrm2] = EPI_distort_fieldmap_image(dwi_slice(:,:,1), -1*deltaB0_slice, echo_space, ...
%                                                  deform_scale_factor, mask_undistorted(:,:,1), mask_dwi_slice(:,:,1));
% 
% dfSize = size(deform_trsfrm1);



% L-B operator 
L = [];
for l = 0:2:shOrder
   L = cat(1, L, l*ones([2*l+1 1]));
end
LB_diag = -1*L.*(L+1);
LB = sparse(diag(LB_diag));
clear LB_diag L


% 1st order finite difference matrix 
finite_diff = createDWithPeriodicBoundary(dwiSize(1), dwiSize(2));
mask_ind = find(mask_undistorted(:,:,1)<=0);
finite_diff(:, mask_ind) = [];
finite_diff([mask_ind; mask_ind+nVoxels],:) = [];

% spatial weighting for spherical smoothing
spherical_weight(mask_ind) = []; % remove unwanted points
spherical_weight = sparse(diag(spherical_weight));

% create reshape index for two PE directions
temp = 1:length(phase_dir);
dir1_index = temp(phase_dir(:)==1);
dir2_index = temp(phase_dir(:)==2);
dir3_index = temp(phase_dir(:)==3);
dir4_index = temp(phase_dir(:)==4);

% % create index for two PE directions
% num_dir_1 = sum(phase_dir(:)==1);
% num_dir_2 = length(phase_dir) - num_dir_1;
% dir1_index = zeros(num_dir_1,1);
% dir2_index = zeros(num_dir_2,1);
% 
% dir1_count = 1;
% dir2_count = 1;
% for iDir = 1:nDirection
%    if phase_dir(iDir) == 1
%       dir1_index(dir1_count) = iDir;
%       dir1_count = dir1_count+1;
%    else
%       dir2_index(dir2_count) = iDir;
%       dir2_count = dir2_count+1;
%    end
% end


% create vector for observed data + zeros for regularization
d_len = ndim3*dfSize(1) + Ysize(2)*dfSize(2) + size(finite_diff, 1)*Ysize(1);
d = zeros(d_len,1);
d(1:sum(mask_dwi_slice(:)>0),:) = dwi_slice(mask_dwi_slice);


Afun = @(x, transp_flag)transform_X_invLeastSq_spharm(x, transp_flag, ...
                                    deform_trsfrm1, deform_trsfrm2, deform_trsfrm3, deform_trsfrm4,...
                                    Y_sh, alpha(1)*LB, alpha(2)*finite_diff, isfull, num_diffusion_dir,...
                                    dir1_index, dir2_index, dir3_index, dir4_index, spherical_weight);

[c, flag, relres]  = lsqr(Afun, d, 1e-6, 1000);

c = reshape(c, dfSize(2), []);
dwi = c * (Y_sh');


% put back estimated data back to original dimension
c_out = zeros([dwiSize(1:2) Ysize(2)]);
data_out = zeros([dwiSize(1:2) Ysize(1)]);

data_out(mask_undistorted(:,:,ones(1, Ysize(1)))) = dwi(:);
c_out(mask_undistorted(:,:,ones(1, Ysize(2)))) = c(:);


end


function y = transform_X_invLeastSq_spharm(x, transp_flag, ...
                         deform_trsfrm1, deform_trsfrm2, deform_trsfrm3, deform_trsfrm4,...
                         Y_sh, LB, finite_diff, isfull, num_diffusion_dir,...
                         dir1_index, dir2_index, dir3_index, dir4_index, spherical_weight)

if strcmp(transp_flag,'transp')  % y = A'*x
   
   LB_length = size(Y_sh, 2)*size(deform_trsfrm1, 2);
   FD_length = size(finite_diff, 1)*size(Y_sh, 1);
   dwi_length = length(x) - (LB_length+FD_length);
   
   x_LB = reshape(x(1+dwi_length:dwi_length+LB_length), size(LB, 1), []);
   y_LB = x_LB' * LB; % LB * (x_LB');
   y_LB = spherical_weight * y_LB; % spherical weight
   
   x_FD = reshape(x(end-FD_length+1:end), 2*size(deform_trsfrm1, 2), []);
   y_FD = (finite_diff' * x_FD) * Y_sh;

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
      
      y_deform_out = (y1+y2+y3+y4)*Y_sh;
   else
      y_deform = [y1 y2 y3 y4];
      clear y1 y2 y3 y4
      
      z = zeros(size(y_deform));
      z(1:end, [dir1_index dir2_index dir3_index dir4_index]) = y_deform;
      clear y_deform
      
      y_deform_out = z * Y_sh;
      clear z
   end
   y = y_deform_out + y_LB + y_FD;
   y = y(:);
   
elseif strcmp(transp_flag,'notransp') % y = A*x
   c = reshape(x, size(deform_trsfrm1, 2), []);
   dwi = c * (Y_sh');
   
   % L-B operator on spharm coeff
   y_LB = LB * (c');
   y_LB = y_LB * spherical_weight; % apply spherical weight

   % B0 deformation
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
      y1 = deform_trsfrm1 * dwi(1:end, dir1_index);
      y2 = deform_trsfrm2 * dwi(1:end, dir2_index);
      y3 = deform_trsfrm3 * dwi(1:end, dir3_index);
      y4 = deform_trsfrm4 * dwi(1:end, dir4_index);
   end
   y_deform = [y1 y2 y3 y4];
   clear y2 y1 y3 y4
   y_deform_out = zeros(size(y_deform));
   y_deform_out(1:end, [dir1_index dir2_index dir3_index dir4_index]) = y_deform;
   clear y_deform
   
   % finite difference
   y_FD = finite_diff * dwi;
   
   y = [y_deform_out(:); y_LB(:); y_FD(:)];
end

end
