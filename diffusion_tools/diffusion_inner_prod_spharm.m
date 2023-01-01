% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [inner_prod, inner_prod_anisotropy] = diffusion_inner_prod_spharm(S1, bMatrices1, S2, bMatrices2, shOrder)
% Computes innner product between two diffusion-weighted images S1 & S2. Last dimension is the
% no of diffusion directions. bMatrices - 3*3*n 
% This function expects to include b=0 image corresponding to zero bmatrices.
%
% inner_prod - inner product using all coefficients
% inner_prod_anisotropy - inner product after discarding first SH coeff (isotropic component)
%

S1Size = size(S1);
S2Size = size(S2);

if ~isequal(S1Size(1:end-1), S2Size(1:end-1))
   error('D1 & D2 should have same no of voxels')
end

S1 = double(S1);
S2 = double(S2);

% throw away negative numbers
S1(S1<0) =  0;
S2(S2<0) =  0;

% Find spharm coeff for D1 & D2
c1 = fit_ADC(S1, bMatrices1, shOrder);
c2 = fit_ADC(S2, bMatrices2, shOrder);

% compute inner products
c1 = reshape(c1, prod(S1Size(1:end-1)), []);
c2 = reshape(c2, prod(S2Size(1:end-1)), []);

% normalize coeffs 
c1 = c1./repmat(sqrt(sum(c1.^2, 2)), [1 size(c1,2)]);
c2 = c2./repmat(sqrt(sum(c2.^2, 2)), [1 size(c2,2)]);


inner_prod = sum(c1.*c2, 2);
inner_prod_anisotropy = sum(c1(:,2:end).*c2(:,2:end), 2);

if numel(S1Size(1:end-1))==1
   S1Size(end) = 1; % add 1 to make it reshape-able
   S1Size(end+1) = 1;
end
inner_prod = reshape(inner_prod, S1Size(1:end-1));
inner_prod_anisotropy = reshape(inner_prod_anisotropy, S1Size(1:end-1));

end

