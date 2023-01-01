% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [ J1 J2 J3 ] = jacobian_nii( data, smoothType )
% Calculates the Jacobian of the affine transform, which can
% be used to calculate transform matrix
%
%   data - 4D matrix, first three dimension denotes location of each voxel.
%          Fourth dimension denoted the x,y,z map location for that voxel.
%          'data' should be in nii format
%
%   smoothType - Defaults to 'smooth'. When set to 'no_smooth' does not smooth
%                the mask before calculating the jacobian.
%
%J1, J2, J3 - Jacobian of each component(x,y,z)
%

[len bre dep dim4] = size(data.img);

% Jacobian without smoothing.
if exist('smoothType', 'var') && strcmpi(smoothType, 'no_smooth')
   J1 = zeros([len bre dep dim4]);
   J2 = zeros([len bre dep dim4]);
   J3 = zeros([len bre dep dim4]);
   
   for m = 1:dim4
      [temp1 temp2 temp3] = gradient(data.img(:,:,:,m));
      J2(:,:,:,m) = temp1;   % order of J1 & j2 are interchanged to accomodate the matlab's convention with ndgrid. check gradient & ndgrid.
      J1(:,:,:,m) = temp2;
      J3(:,:,:,m) = temp3;
      % fprintf('jacobian: %2.2f %%\n',m/dim4*100)
   end
   clear J
   
% Jacobian with smoothing
elseif exist('smoothType', 'var') && strcmpi(smoothType, 'smooth')

   J1 = zeros([len bre dep dim4]);
   J2 = zeros([len bre dep dim4]);
   J3 = zeros([len bre dep dim4]);
   
   for m = 1:dim4
      [temp1 temp2 temp3] = gradient(smooth3(data.img(:,:,:,m),'gaussian',[7 7 7],2));
      J2(:,:,:,m) = temp1;   % order of J1 & j2 are interchanged to accomodate the matlab's convention with ndgrid. check gradient & ndgrid.
      J1(:,:,:,m) = temp2;
      J3(:,:,:,m) = temp3;
      % fprintf('jacobian: %2.2f %%\n',m/dim4*100)
   end
   
else
   error('Specify the smoothType parameter correctly.');
end

end

