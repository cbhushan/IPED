% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


% This function returns the SHORE basis.
function [ basis,R,Snew,N ] = shore_basis(qrad, qcart,rad_ord,zeta) 
    basis = [];
    % Angular part of the basis - Spherical harmonics
    [S,L] = sph_harm_basis(qcart,rad_ord,2);
    Snew=[];
    ind = 1;
    for n = 0:rad_ord, %Radial order
        nInd = 1;
        for l = 0:2:n, % Ang order
            Snew =[Snew S(:,L==l)];
            for m = -l:l,
                % Radial part
                R(:,ind) = ((2*factorial(n-l))/(zeta^1.5*gamma(n+1.5)))^(0.5)*(qrad.^2/zeta).^(l/2).*exp(-(qrad.^2)/(2*zeta)).*...
                    laguerrePoly(n-l,l+0.5,(qrad.^2/zeta));
                N(ind) = nInd;
                ind = ind+1;
                nInd = nInd+1;
            end;
        end;
    end;
    basis = R.*Snew;    
end
