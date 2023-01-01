% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [c, Y] = fit_ADC(dwi, bMatrices, shOrder)
% Estimates coefficients of shOrder of modified SH. 
% 
% dwi - Diffusion weighted data, with b=0 image(s). Last dimension is the no of diffusion
%     directions(n) . eg: X*Y*Z*n
%
% bMatrix - 3*3*n
%
% shOrder - Order of SH basis
%
% c - coefficients of first nTerms of modified SH. Last dimension is the
%     estimated coefficients. Eg: X*Y*Z*nTerms
% Y - First nTerms modified SH basis
%

dwiSize = size(dwi);
nVoxels = prod(dwiSize(1:end-1));
nDirection = size(bMatrices, 3);

dwi = reshape(dwi, nVoxels, [])';

if nDirection~=size(dwi,1)
   error('Number of directions in dwi and bMatrices does not match.');
end

if sum(dwi(:)<0) > 0
   warning('Negative values are found in D. Result may be complex valued.')
end

% find b=0 images or zero bmatrices 
temp = reshape(bMatrices, 9, []);
temp = sum(temp.^2, 1);
b0mask = (temp==0);
clear temp

% extract b=0 data
if ~any(b0mask)
   error('No zero matrices (b=0 image) found in D')
   
elseif sum(b0mask)==1
   s0 = dwi(b0mask, :);
   
elseif sum(b0mask)>1
   s0 = mean(dwi(b0mask, :), 1);
   disp([num2str(sum(b0mask)) ' zero b-matrice(s) are found. Using mean of all b=0 image as b=0 image.'])
end
dwi(b0mask, :) = [];
bMatrices(:,:,b0mask) = [];
clear b0mask


% get eigen val/vec and normalize D
nDirection = size(bMatrices, 3);
vec = zeros(nDirection, 3);
adc = zeros(size(dwi));
for iDir = 1:nDirection
   [V, bval] = sorted_eig(bMatrices(:,:,iDir));
   vec(iDir,:) = V(:,1);
   adc(iDir,:) = -1 * log(dwi(iDir,:)./s0) / bval(1,1);
end
clear s0


% modified SH basis
[theta, phi, ~] = cart2sph(vec(:,1),vec(:,2),vec(:,3));
[Y, nTerms] = spharm_Hardi(shOrder, theta, phi); % basis matrix

c = Y \ adc;
c = reshape(c', [dwiSize(1:end-1) nTerms]);

end

