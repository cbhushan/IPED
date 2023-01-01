% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [ v ] = cross_vec( a, b )
% Computes cross product axb, where vectors are stored in 4th
% dimesion of matrix a & b (usually the eigenvectors)
%
% a, b - 4D matrix (n x m x k x 3)
%
% Should be faster than matlab implementation in some cases.

v = zeros(size(a));

v(:,:,:,1) = (a(:,:,:,2).*b(:,:,:,3)) - (a(:,:,:,3).*b(:,:,:,2));
v(:,:,:,2) = (a(:,:,:,3).*b(:,:,:,1)) - (a(:,:,:,1).*b(:,:,:,3));
v(:,:,:,3) = (a(:,:,:,1).*b(:,:,:,2)) - (a(:,:,:,2).*b(:,:,:,1));

end

