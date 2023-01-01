% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function Y = spharm(shDegree, shOrder, theta, phi)
% This function generates the Spherical Harmonics basis functions of degree
% 'shDegree' and order 'shOrder' at discrete polar coordinates 
% (thetaRadians, phiRadians).
%
% thetaRadians, phiRadians should have same dimensions and follow matlab
% convention. See sph2cart. 
%
%    phi: elevation angle, interval [-pi/2,pi/2]
%    theta: azimuth angle, interval [0, 2*pi]
%
%  Y - sampled spherical harmonics on the sphere. Same dimension as
%      thetaRadians and phiRadians. 
%


if min(phi(:))<-pi/2 || max(phi(:))>pi/2
   error('Check range of phiRadians');
end

if shDegree<abs(shOrder)
   error('The ORDER (shOrder) must be less than or equal to the DEGREE(shDegree).');
end

if isvector(phi)
   phi = phi(:)';
   theta = theta(:)';
end

P = legendre(shDegree, sin(phi));

if shDegree ~= 0
   P = squeeze(P(abs(shOrder)+1,:,:));
   if shOrder < 0
      P = P * (-1)^(shOrder) * factorial(shDegree+shOrder)/factorial(shDegree-shOrder);
   end
end

Y = sqrt((2*shDegree+1)*factorial(shDegree-shOrder)/((4*pi)*factorial(shDegree+shOrder)))...
    * (P.*exp(1i*shOrder*theta));
 
end
