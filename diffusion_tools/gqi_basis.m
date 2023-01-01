% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


% This function return the basis for GQI reconstruction
function [ basis ] = gqi_basis(sigma,l_delta,bval,qcart,pts,type,del_t)
if(size(qcart,1) == size(bval,1))
    sx = 1;
else
    sx = size(qcart,1);
end;

Aq = 4*pi*bval/(4*pi^2*del_t);
ld_ang = repmat(sqrt(bval*0.018),[sx size(pts',2)]); %0.01506
%ld_ang = repmat(sqrt(bval*0.018),[sx 3]); %0.01506
%qcart = qcart.*ld_ang;

if(type == 1)
     basis = repmat(Aq,[1 size(pts',2)])*l_delta.*sinc(sigma*ld_ang.*(qcart*pts')/pi);
 %basis = sinc(sigma*(qcart*pts')/pi);
else
    tic
    x=sigma*ld_ang.*(qcart*pts');
    basis = (2*x.*cos(x) + (x .* x - 2).*sin(x))./(x.^3); %*x.*x);
    basis(x == 0) = 1/3;
    basis = (l_delta.^3)*repmat(Aq,[1 size(basis,2)]).*basis;
    
end;
end

