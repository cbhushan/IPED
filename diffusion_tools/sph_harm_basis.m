% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [S,L,Z] = sph_harm_basis(SamplingLocations,HarmonicOrder,type)
% Justin Haldar (jhaldar@usc.edu) 09/02/2011
%
% SamplingLocations -- an Nx3 matrix of vectors describing the
%                      locations where the spherical harmonics are to be sampled
% HarmonicOrder -- the maximum order of the spherical harmonic representation
% type -- 1:  the standard complex spherical harmonics
%         2:  the real symmetric spherical harmonics
%
% S -- an NxM matrix of sampled spherical harmonics, where M is the number
%      of basis functions 
% L -- Mx1 vector describing the order of each spherical harmonic basis
%      function
% Z -- Mx1 vector describing the other harmonic

N = size(SamplingLocations,1);

%% Standard Complex Spherical Harmonics

[theta,phi,rho] = cart2sph(SamplingLocations(:,1),SamplingLocations(:,2),SamplingLocations(:,3));
phi = -phi+pi/2;

if (type ==1)
    S = zeros(N,HarmonicOrder^2+2*HarmonicOrder+1);
    S(:,1) = 1/sqrt(4*pi);
        L = zeros(HarmonicOrder^2+2*HarmonicOrder+1,1);
    for ell=1:HarmonicOrder
        Pell=legendre(ell,cos(phi))';
        Pell = (ones(N,1)*(sqrt((2*ell+1)/(4*pi)*factorial(ell-[-ell:ell])./factorial(ell+[-ell:ell])))).*[flipdim(Pell(:,2:end)*diag((-1).^[1:ell].*factorial(ell-[1:ell])./factorial(ell+[1:ell])),2),Pell].*exp(complex(0,1)*theta*[-ell:ell]);
        S(:,1 + ell^2:1+ell^2+2*ell) = Pell;
            L(1 + ell^2:1+ell^2+2*ell) = ell;
    end
elseif (type==2)
    S = zeros(N,1/2*(HarmonicOrder+1)*(HarmonicOrder+2));
    S(:,1) = 1/sqrt(4*pi);
    index = 1;
    L = zeros(1/2*(HarmonicOrder+1)*(HarmonicOrder+2),1);
    Z = zeros(1/2*(HarmonicOrder+1)*(HarmonicOrder+2),1);
    for ell = 2:2:HarmonicOrder
        Pell=(ones(N,1)*(sqrt((2*ell+1)/(4*pi)*factorial(ell-[0:ell])./factorial(ell+[0:ell])))).*legendre(ell,cos(phi))'.*exp(complex(0,1)*theta*[0:ell]);
        Pell = [flipdim(sqrt(2)*real(Pell(:,2:end)),2),Pell(:,1),sqrt(2)*(ones(N,1)*(-1).^[2:ell+1]).*imag(Pell(:,2:end))];
        S(:,1+ell*(ell-1)/2:(ell+1)*(ell+2)/2) = Pell;
            L(1+ell*(ell-1)/2:(ell+1)*(ell+2)/2) = ell;
            Z(1+ell*(ell-1)/2:(ell+1)*(ell+2)/2) = [-ell:ell];
    end
    
end

