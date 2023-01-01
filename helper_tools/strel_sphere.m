% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function se = strel_sphere(sz)
[x1,y1,z1] = ndgrid(-ceil(sz):ceil(sz));
se = (sqrt(x1.^2 + y1.^2 + z1.^2) <=sz);
end
