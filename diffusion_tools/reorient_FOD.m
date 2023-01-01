% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [cTransform] = reorient_FOD(c, F)
% Estimates coefficients of shOrder of modified SH. 
% 
% c - coefficients of modified SH. Last dimension is the
%     estimated coefficients. A 4D Eg: X x Y x Z x nTerms
%
% bMatrix - 3 x 3 x nDir . Last dimension is the number of dimensions
%
% F - Local affine transform matrix. 3 x 3 (applies same transform to all voxels)

nDirections = 300;

sizeC = size(c);
nTerms = sizeC(end);

if length(sizeC) > 1
   nVoxels = prod(sizeC(1:end-1));
else
   nVoxels = 1;
end

c = reshape(c, nVoxels, [])';
shOrder = (sqrt(1+8*nTerms)-3)/2;

load('300dir.mat');
theta = theta';
phi = phi';
%[theta, phi] = generate_uniform_dir(nDirections);

M = spharm_Hardi(shOrder, theta, phi);
w = M \ c; % estimate weights

% Apply transform
[vecX, vecY, vecZ] = sph2cart(theta, phi, 1);
vec = [vecX vecY vecZ];
vecT = vec * F';

[thetaT, phiT, ~] = cart2sph(vecT(:,1),vecT(:,2),vecT(:,3));
MTransform = spharm_Hardi(shOrder, thetaT, phiT);

cTransform = MTransform * w;
cTransform = reshape(cTransform', sizeC);

end
