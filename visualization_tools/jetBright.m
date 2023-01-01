% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [c, cntr] = jetBright(m)
% A colormap which chops-off part of the top & bottom of jet() colormap where it starts getting
% darker. Works well with m>26. Does NOT returns colormap with size m - use with caution!

if nargin < 1
   m = size(get(gcf,'colormap'),1);
end

cjet = jet(m);

% this sum should be 1 in all area except top and bottom end of jet()
tt = cjet(:,1) + cjet(:,3);

cntr = 0;
while cntr<(m/3) && max(tt(cntr+1), tt(end-cntr))<0.98
   cntr = cntr + 1;
end

% maintain symmetry
c = cjet(1+cntr:end-cntr, :);

end
