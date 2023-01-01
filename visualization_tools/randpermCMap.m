% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function c = randpermCMap(cmap, type)
% Randomly permute colormap and adds a black element at the lower-end.
sz = size(cmap);
if nargin==1 || type == 1
   p = randperm(sz(1));   
elseif type == 2
   n = ceil(sz(1)/2);
   t1 = [1:n; n+1:2*n];
   p = t1(1:sz(1));
end
c = cat(1, [0.2 0.2 0.2], cmap(p,:));
end

