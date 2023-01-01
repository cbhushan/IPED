% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


% This function returns the fourier transform of the SHORE basis. This is
% the EAP basis corresponding to SHORE basis
function [ basis,R,Snew,R1,R2] = shoredual_basis(qrad, qcart,rad_ord, ang_ord,zeta)   
     basis = [];
    % Angular part of the basis - Spherical harmonics
    [S,L] = sph_harm_basis(qcart,rad_ord,2);
    Snew=[];
    ind = 1;
    for n = 0:rad_ord, %Radial order
        for l = 0:2:n, % Ang order
            Snew =[Snew S(:,L==l)];
            for m = -l:l,
                % Radial part
                p = n - l;
                c=(4*pi^2*zeta);
                R(:,ind) = ((-1)^(n - l/2))*sqrt((2*(c^(1.5))*factorial(p))/(gamma(n+1.5)))*((c*qrad.^2).^(l/2)).*...
                    exp(-((c/2)*(qrad.^2))).*laguerrePoly(p,l+0.5,(c*qrad.^2));
                R2(:,ind) = ((-1)^(n - l/2))*sqrt((2*(c^(1.5))*factorial(p))/(gamma(n+1.5)))*((c*qrad.^2).^(l/2)).*...
                    exp(-((c/2)*(qrad.^2)));
                R1(:,ind) = laguerrePoly(p,l+0.5,(c*qrad.^2));
                ind=ind+1;
            end;
        end;
    end;
    basis = R.*Snew;
end
