% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function d = extrapolateInVolSmooth(vol, knw_msk, res, lmbd, xinit, tol, maxit)
% Extrapolates data outside the knw_msk in the 3D volume maintaining the smoothness. 
% Values inside the knw_msk may NOT be preserved. Also see extrapolateInVol.m 
% res is the resolution/voxel size along each dimension. 
% 
% Can be VERY memory intensive!

sz = size(vol);
if ~isequal(sz, size(knw_msk))
   error('size of vol and mask does not match')
end
vol = double(vol);

if ~exist('tol', 'var')
   tol = 1e-7;
end

if ~exist('maxit', 'var')
   maxit = 1000;
end

if ~exist('xinit', 'var')
   xinit = ones(prod(sz), 1);
end

knw_msk = knw_msk>0;

mask_knw = createMaskOperators(knw_msk);
D = createDNoBoundary3DresWt(sz(1), sz(2), sz(3), 1./res);
b = [mask_knw*vol(:); zeros(size(D,1),1)];
A = [mask_knw; lmbd*D];
clear D mask_knw
x = lsqr(A, b, tol, maxit, [], [], xinit(:));
d = reshape(x, sz);

end

