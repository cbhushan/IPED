% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function c = pinkLowHigh(m)
% A colormap which is similar to pink on both low and high end of color scale. 

if nargin < 1
   m = size(get(gcf,'colormap'),1); 
end

cpink = pink(round(m/2));
c = [cpink; flipdim(cpink(:,[2 1 3]) ,1)];

end
