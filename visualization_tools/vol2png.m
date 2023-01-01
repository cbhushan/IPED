% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function vol2png(vol, clim, pngFileBase, type)
% Writes a png image (montage) of 3D volume.
%   vol - 3D/4D volume. See type
%   clim - min-max range of intensities
%   pngFileBase - Output filebase with full path. Saved files are named as [pngFileBase '.1.png']
%   type - Can be either 'rgb' or a colormap. When set to 'rgb', vol must be a 4D volume (x*y*z*3).
%          To use a colormap, type must be a Nx3 matrix of double representing a RGB colormap. 
%          Ex: 
%             vol2png(Imv, [0 1], 'out', jet); % It will write out png files with jet colormap. 
%

if exist('type', 'var')
   if ischar(type) && ~strcmpi(type, 'rgb')
      error('type can be either ''rgb'' or Nx3 matrix for colormap')
   end
   
   if isnumeric(type) && size(type,2) ~= 3
      error('type can be either ''rgb'' or Nx3 matrix for colormap')
   end   
else
   type = [];
end

if ischar(type) && size(vol, 4)~=3
   error('vol must be 4D volume (x*y*z*3) for type RGB.')
elseif ~ischar(type) && ndims(vol)~=3
   error('vol must be 3D volume (x*y*z)')
end
vol = double(vol);

slices3 = 1:size(vol, 3);
slices1 = 1:size(vol, 1);
slices2 = 1:size(vol, 2);
nRows = NaN;
nCols = NaN;

img = (vol-clim(1))/(clim(2)-clim(1));
clim = [0 1];
img(img>1) = 1;
img(img<0) = 0;
clear vol

% Apply colormap
if ~isempty(type) && ~ischar(type)
   img_c = zeros([size(img) 3]);
   crange = linspace(0, 1, size(type, 1)); 
   for k = 1:3
      img_c(:,:,:,k) = interp1(crange, type(:,k), img);
   end
   img = img_c; 
   clear img_c
end

E = permute(img, [2 1 4 3]);
E = flipdim(E, 1);
I1 = montage_image(E, 'size', [nRows nCols], 'Indices', slices3, 'DisplayRange', clim);
imwrite(I1, [pngFileBase '.1.png'], 'png');
clear I1

E = permute(img, [3 1 4 2]);
E = flipdim(E, 1);
I2 = montage_image(E, 'size', [nRows nCols], 'Indices', slices2, 'DisplayRange', clim);
imwrite(I2, [pngFileBase '.2.png'], 'png');
clear I2

E = permute(img, [3 2 4 1]);
E = flipdim(E, 1);
E = flipdim(E, 2);
I3 = montage_image(E, 'size', [nRows nCols], 'Indices', slices1, 'DisplayRange', clim);
imwrite(I3, [pngFileBase '.3.png'], 'png');
clear I3 E A img

end

