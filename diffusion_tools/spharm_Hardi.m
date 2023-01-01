% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [ Y, nTerms, L, Z ] = spharm_Hardi(shOrder, theta, phi)
% Generates modified spherical harmonic basis for HARDI signal (real and symmetric) truncated
% after order 'shOrder'. 
%
% thetaRadians, phiRadians should have same dimension. They should follow
% matlab convention. See sph2cart
%
%    length(thetaRadians) = length(phiRadians) = no. of direction of HARDI
%    acquisition.
%
%    phi: elevation angle, interval [-pi/2,pi/2]
%    theta: azimuth angle, interval [0, 2*pi]
%
%  Y - Sampled spherical harmonics on sphere. Last dimension is of length nTerms. If the input
%      is a vector then output is always a 2D matrix.  
%
%  nTerms - number of spharm basis functions
%  L - vector of length nTerms describing the order of each spherical harmonic basis function
%  Z - vector of length nTerms describing the degree of each spherical harmonic basis function
%

if min(phi(:))<-pi/2 || max(phi(:))>pi/2
   error('Check range of phiRadians');
end

thetaSize = size(theta);
if length(thetaSize)==2 && min(thetaSize)==1
   thetaSize = max(thetaSize);
end

shOrder = floor(shOrder/2)*2; % greatest even integer <= shOrder
nTerms = (shOrder+1)*(shOrder+2)/2;

theta = theta(:);
phi = phi(:);

N = length(phi);
Y = zeros([N, nTerms]);
L = zeros(nTerms,1);
Z = zeros(nTerms,1);

% faster (vectorized) implementation copied from Justin
Y(:,1) = 1/sqrt(4*pi);
for ell = 2:2:shOrder
   Pell=(ones(N,1)*(sqrt((2*ell+1)/(4*pi)*factorial(ell-(0:ell))./factorial(ell+(0:ell))))) .* legendre(ell,sin(phi))' .* exp(1i*theta*(0:ell));
   Pell = [flipdim(sqrt(2)*real(Pell(:,2:end)),2),  Pell(:,1),  sqrt(2)*(ones(N,1)*(-1).^(2:ell+1)).*imag(Pell(:,2:end))];
   Y(:,1+ell*(ell-1)/2:(ell+1)*(ell+2)/2) = Pell;
   L(1+ell*(ell-1)/2:(ell+1)*(ell+2)/2) = ell;
   Z(1+ell*(ell-1)/2:(ell+1)*(ell+2)/2) = (-ell:ell);
end

% Old slower implementation 
% for l = 0:2:shOrder
%    for m = -l:l
%       if m<0 
%          Y(:, iTerm) = sqrt(2)*real(spharm(l, abs(m), theta, phi));
%       elseif m>0
%          Y(:, iTerm) = sqrt(2)*(-1)^(m+1)*imag(spharm(l, m, theta, phi));
%       else
%          Y(:, iTerm) = spharm(l, m, theta, phi);
%       end
%      iTerm = iTerm+1;
%    end
% end

Y = reshape(Y, [thetaSize, nTerms]);

end
