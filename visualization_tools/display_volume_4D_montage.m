% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function I = display_volume_4D_montage(Img, varargin)
% montage display of 4D volumes. 
%  usage:
%      display_volume_4D_montage(Img)
%      display_volume_4D_montage(Img, dim)           % default dim=4
%      display_volume_4D_montage(Img, clim)          % clim = [low high]
%      I = display_volume_4D_montage(Img, ...)       % Returns montaged 3D volume
%

if ndims(Img)~=4
   error('Image must be 4D volume')
end

for n = 2:nargin
   if length(varargin{n-1}) == 2
      clim = varargin{n-1};
   else
      dim = varargin{n-1};
   end
end

if ~exist('clim', 'var')
   I_s = sort(Img(:), 'ascend');
   low = I_s(max(floor(length(I_s)*0.02), 1));
   high = I_s(floor(length(I_s)*0.985));
   clim = [low high];
   
   if clim(1)==clim(2)
      clim = clim + [-1 1];
   end
end

if ~exist('dim', 'var')
   dim = 4;
end

switch dim
   case 4
      Img = permute(Img, [1 2 4 3]);
   case 3
      % do nothing
   case 2
      Img = permute(Img, [1 3 2 4]);
   case 1
      Img = permute(Img, [3 2 1 4]);
   otherwise
      error('dim can not be more than 4');
end

sz = size(Img);
nRows = ceil(4*sqrt(sz(4)/12));
nCols = ceil(3*sqrt(sz(4)/12));
Img = flipdim(Img, 1);
I = zeros(sz(1)*nRows, sz(2)*nCols, sz(3));
for k = 1:sz(3)
   I(:,:,k) = montage_image(Img(:,:,k,:), 'size', [nRows nCols]);
end
clear Img

if nargout<1  
   display_volume(I, clim);
   clear I;
end

end

