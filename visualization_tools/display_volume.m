% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function h = display_volume(Img, varargin)
% Wrapper arround imdisp to display a 3D volume in intractive figure
% window. In case of 4D input it must be a RGB 3D volume.
%  usage:
%      display_volume(Img)
%      display_volume(Img, dim)           % default dim=3
%      display_volume(Img, clim)          % clim = [low high]
%      display_volume(Img, type)          % type can be 'm' or 'rgb'
%      h = display_volume(Img, ...)
%

Img = double(Img);

mont = false;
rgb = false;

if ndims(Img)==4
   if size(Img, 4)==3
      rgb = true;
   else
      error('4D volumes must be RGB volumes with x*y*z*3.')
   end
end

for n = 2:nargin
   if ischar(varargin{n-1}) && strcmpi(varargin{n-1}, 'm')
      mont = true;
   elseif ischar(varargin{n-1}) && strcmpi(varargin{n-1}, 'rgb')
      rgb = true;
   elseif isnumeric(varargin{n-1}) && length(varargin{n-1}) == 2
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
   dim = 3;
end

clim
if rgb
   Img = (Img-clim(1))/(clim(2)-clim(1));
   clim = [0 1];
end

if rgb && ndims(Img)==3
   Img = permute(Img(:,:,:,[1 1]), [1 2 4 3]);
end

switch dim
   case 3
      E = permute(Img, [2 1 4 3]);
   case 2
      E = permute(Img, [3 1 4 2]);
   case 1
      E = permute(Img, [3 2 4 1]);
   otherwise
      error('dim can not be more than 3');
end

E = flipdim(E, 1);
h = figure;

if mont
   %imdisp(E, gray(256), clim);
   montage(E, 'DisplayRange', clim);
else
   imdisp(E, gray(256), clim, 'Size', 1);
end

end

