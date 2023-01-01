% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function d = extrapolateInVol(vol, knw_msk, res, boundary_sz, xinit, tol, maxit)
% Extrapolates data outside the knw_msk in the 3D volume. Values inside the knw_msk would be
% exactly preserved (values at the boundary of mask is mainly used for extrapolation.) 
% res is the resolution/voxel size along each dimension. 
% boundary_sz is the size boundary to use from the know data (typically set to 2).
% 
% Can be VERY memory intensive!

sz = size(vol);
if ~isequal(sz, size(knw_msk))
   error('size of vol and mask does not match')
end

if ~exist('tol', 'var')
   tol = 1e-7;
end

if ~exist('maxit', 'var')
   maxit = 1000;
end

knw_msk = knw_msk>0;
vol = double(vol);

[x1,y1,z1] = ndgrid(-boundary_sz:boundary_sz);
se = (sqrt(x1.^2 + y1.^2 + z1.^2) <=boundary_sz);

[mask_knw, unmask_knw, ~, unmask_unknw] = createMaskOperators(knw_msk);
D = createDwithMask(sz(1), sz(2), sz(3), imdilate(~knw_msk, se), 1./res);
b = -1*D*unmask_knw*mask_knw*vol(:);
A = D*unmask_unknw;
clear D

if ~exist('xinit', 'var') || isempty(xinit)
   xinit = ones(size(unmask_unknw,2), 1);
end

x = lsqr(A, b, tol, maxit, [], [], xinit);

d = unmask_knw*mask_knw*double(vol(:)) + unmask_unknw*x;
d = reshape(d, sz);

end

