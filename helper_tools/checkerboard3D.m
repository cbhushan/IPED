% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function I = checkerboard3D(n, varargin)
% Creates a 3D checkerboard. 
% Usage: 
%    I = checkerboard3D(n) % nxnxn size
%    I = checkerboard3D([n1 n2 n3]) % n1xn2xn3 size
%    I = checkerboard3D(n1, n2, n3) % n1xn2xn3 size
%

if nargin==1 && length(n) == 1
   out_sz = [n n n];
elseif nargin==1 && length(n) == 3
   out_sz = n;
elseif nargin==3
   out_sz = [n varargin{1} varargin{2}];
else
   error('Incorrect input. See usage.')
end

sz = ceil(out_sz/2)*2; % even

v1 = zeros(sz(1:2));
v1(1:2:end, 1:2:end) = 1;
v1(2:2:end, 2:2:end) = 1;

v2 = circshift(v1, [1 0]);

v = cat(3,v1,v2);
I = repmat(v, [1 1 sz(3)/2]);
I = I(1:out_sz(1), 1:out_sz(2), 1:out_sz(3));

end
