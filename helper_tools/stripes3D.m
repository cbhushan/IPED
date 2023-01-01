% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function I = stripes3D(n)
% Creates a 3D checkerboard. 
% Usage: 
%    I = stripes3D(n) % nxnxn size
%    I = stripes3D([n1 n2 n3]) % n1xn2xn3 size


if nargin==1 && length(n) == 1
   sz = [n n n];
elseif nargin==1 && length(n) == 3
   sz = n;
else
   error('Incorrect input. See usage.')
end

v1 = zeros(sz(1:2));
v1(:, 1:2:end) = 1;

v2 = zeros(sz(1:2));
v2(1:2:end, :) = 1;

v = cat(3,v1,v2);
I = repmat(v, [1 1 ceil(sz(3)/2)]);
I = I(1:sz(1), 1:sz(2), 1:sz(3));

end
