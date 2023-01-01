% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [cTransform] = reorient_FOD2(c, F)

nDirections = 300;

cSize = size(c);
nVoxels = prod(cSize(1:end-1));
nTerms = cSize(end);
c = reshape(c, nVoxels, [])';
shOrder = (sqrt(1+8*nTerms)-3)/2;

load('300dir.mat');
theta = theta';
phi = phi';

Y = squeeze(spharm_Hardi(shOrder, theta, phi));
odf = Y' * c;

% Apply transform
[vecX, vecY, vecZ] = sph2cart(theta, phi, 1);
vec = [vecX vecY vecZ];
vecT = vec * F';
[thetaT, phiT, ~] = cart2sph(vecT(:,1),vecT(:,2),vecT(:,3));

Y = spharm_Hardi(shOrder, thetaT, phiT); % basis matrix

cTransform = Y' \ odf;
cTransform = reshape(cTransform', cSize);

end
