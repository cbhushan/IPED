% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [c, Y] = fitDWIs(dwi, bMatrices, shOrder, LB_alpha)
% Fits DWIs using SH basis functions. Should be single-shell diffusion
% data - No checks!
%
% dwi - dxn 2D matrix; d is number of diffusion weightings, n is number of voxels
% bMatrices - 3*3*d
% shOrder - Order of SH basis
% LB_alpha - (optional) regularization parameter for L-B operator 
%
% c - nTerms*n matrix. coefficients of modified SH
% Y - First nTerms modified SH basis
%

[nDirection, nVoxels] = size(dwi);

if nDirection~=size(bMatrices,3)
   error('Number of directions in dwi and bMatrices does not match.');
end

% throw away b=0 images
DEout = checkDiffusionEncodingScheme(bMatrices);
vec = zeros(nDirection, 3);
for iDir = 1:nDirection
   [V, ~] = sorted_eig(bMatrices(:,:,iDir));
   vec(iDir,:) = V(:,1);
end
vec(DEout.zero_bval_mask, :) = [];
dwi(DEout.zero_bval_mask, :) = [];
nDirection = sum(~DEout.zero_bval_mask);

% modified SH basis
[theta, phi, ~] = cart2sph(vec(:,1),vec(:,2),vec(:,3));
[Y, nTerms, L] = spharm_Hardi(shOrder, theta, phi); % basis matrix

% fit SH
if exist('LB_alpha', 'var') && LB_alpha>0 
   LB_diag = -1*L.*(L+1);
   LB = LB_alpha * sparse(diag(LB_diag));
   clear LB_diag L
   
   % create vector for observed data + zeros for regularization
   d_len = (nDirection + nTerms) * nVoxels;
   d = zeros(d_len,1);
   d(1:nDirection*nVoxels, :) = dwi(:);
   clear dwi
   
   Afun = @(x, transp_flag)transform_X(x, transp_flag, Y, LB, nDirection, nVoxels, nTerms);   
   c = lsqr(Afun, d, 1e-6, 1000);
   c = reshape(c, nTerms, []);
   
else % no regularization
   c = Y \ dwi;   
end

end


function y = transform_X(x, transp_flag, Y, LB, nDirection, nVoxels, nTerms)

if strcmp(transp_flag,'transp')  % y = A'*x
   n = nDirection*nVoxels;
   dwi = reshape(x(1:n), nDirection, []);
   y1 = (Y') * dwi;
   
   % L-B operator
   c_res = reshape(x(1+n:end), nTerms, []);
   y2 = LB * c_res; % LB is diagonal matrix 
   
   y = y1+y2;
   y = y(:);
   
elseif strcmp(transp_flag,'notransp')  % y = A*x
   c = reshape(x, nTerms, []);
   dwi = Y * c;
   y2 = LB * c;
   
   y = [dwi(:); y2(:);];
end

end



