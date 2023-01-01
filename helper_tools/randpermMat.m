% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function B = randpermMat(A)
% Randomly permute matrix A along 1st dim

if numel(A)==1
   B = A;
   return
end

sz = size(A);
p = randperm(sz(1));
A = reshape(A, sz(1), []);
B(p, :) = A;
B = reshape(B, sz);

end
