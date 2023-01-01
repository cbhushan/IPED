% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [R, N] = shore_odf_factor(zeta,rad_ord)
    ind = 1;
    cz = 2*(4*pi^2*zeta)^1.5;
    for n = 0:rad_ord, %Radial order
        nInd = 1;
        for l = 0:2:n, % Ang order
            for m = -l:l,
                % The analytical formula of the integral : Eq C.2, Merlet,
                % 2013
                R(:,ind) = (((-1)^(n - l/2))/cz)*(((cz*factorial(n-l))/gamma(n+1.5))^0.5)*((gamma(l/2+1.5)*gamma(n+1.5))/(gamma(l+1.5)*factorial(n-l)))*(0.5^(-l/2 - 1.5))*(HyperGeometric2F1([-n+l,l/2 + 1.5],l+1.5,2));
                  %R(:,ind) = (((-1)^(n - l/2))/cz)*(((cz*factorial(n-l))/gamma(n+1.5))^0.5)*((gamma(l/2+1.5)*gamma(n+1.5))/(gamma(l+1.5)*factorial(n-l)))*(0.5^(-l/2 - 1.5))*(taylora2f1(-n+l,0,l/2 + 1.5,0,l+1.5,0,2,0,1e-6));
                N(ind) = nInd;
                ind = ind+1;
                nInd = nInd+1;
            end;
        end;
    end;
%     basis = repmat(R,[size(Snew,1),1]).*Snew; 
end
