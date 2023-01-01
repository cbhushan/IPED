% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function cout = randCMapHue(sz, type, n_steps)
% Tries to generate high hue-contrast colormap. This is basically a smarter & better version of hsv
% colormap with few more options.
%  Usage:
%     cout = randCMapHue(sz)
%     cout = randCMapHue(sz, type)
%     cout = randCMapHue(sz, type, n_steps)
%
%   sz - scalar size of desired colormap.
%   type - 'randperm' (default), 'adjperm' or 'noperm'. 'adjperm' tries to make adjacent colors to
%           be high contrast colors (like red -> green -> blue). 'noperm' applies no permutation to
%           colors.
%   n_steps - Number of saturation levels away from full color saturation. Default value is based
%             on sz. This changes the way color-space is transversed - see the code below. 
%             (Use value of 1 to get the behavior of previous version of this function)
%
% Note that randomness is this function is limited to shuffling of selected high-contrast hue. The
% selected set of colors will be same for a fixed sz. Further, it is difficult to generate
% high-contrast colormap with lots of color (eg. sz>50) - try and experiment with n_steps. 
%
% % Following code plots full HSV colorspace
% [h, s, v] = ndgrid(linspace(0,1,100), linspace(0,1,100), linspace(0,1,20));
% sz = size(h);
% rgbimg = reshape(hsv2rgb([h(:) s(:) v(:)]), [sz 3]);
% display_volume(rgbimg, 3)
% display_volume(rgbimg, 1)

if ~exist('type', 'var')
   type = 'randperm';
end

if ~exist('n_steps', 'var')
   if sz<=10
      n_steps = 1;
   elseif sz<26
      n_steps = 2;
   elseif sz<50
      n_steps = 3;
   else
      n_steps = 4;
   end
end

% factor by which adjacent colors should be adjusted
% (really useful in case of lots of colors)
if n_steps < 1
   s_scl = ones(1, sz);
   v_scl = ones(1, sz);
else
   v_min = 0.75;
   s_min = 0.7;
   v_scl = linspace(v_min, 1, n_steps+1);
   v_scl = repmat(flipdim(v_scl, 2), [1 ceil(sz/(n_steps+1))]);
   v_scl = v_scl(1:sz);
   s_scl = linspace(s_min, 1, n_steps+1);
   s_scl = repmat(flipdim(s_scl, 2), [1 ceil(sz/(n_steps+1))]);
   s_scl = s_scl(1:sz);
end

h = linspace(0, 1, (sz+1));
h(end) = []; % throw away hue 1, which is same as hue 0

v = v_scl;
s = s_scl;

c = hsv2rgb([h(:) s(:) v(:)]);

if strcmpi(type, 'randperm')
   p = randperm(sz);
   
elseif strcmpi(type, 'adjperm')
   ind = [1:3:sz 2:3:sz 3:3:sz];
   p = zeros(sz,1);
   p(ind) = 1:sz;
   
elseif strcmpi(type, 'noperm')
   p = 1:sz;
   
else
   error('Unknown type: %s', type)
end

cout = c(p,:);

end
