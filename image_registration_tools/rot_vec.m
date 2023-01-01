% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [ v_rot ] = rot_vec( v, axis, cosTheta )
% Rotates vector v about axis by angle theta. Function takes cosine of
% theta as argument. axis may not be a unit vector. Here vectors are stored
% in 4th dimesion of matrix a & b (usually the eigenvectors) 
% 
%  Uses Rodrigues' rotation formula

axis = norm_vec(axis);

term1 = zeros(size(v));
term1(:,:,:,1) = v(:,:,:,1).*cosTheta;
term1(:,:,:,2) = v(:,:,:,2).*cosTheta;
term1(:,:,:,3) = v(:,:,:,3).*cosTheta;


term2 = zeros(size(v));
sinTheta = sqrt(1 - (cosTheta.^2));
temp = cross_vec(axis, v);
term2(:,:,:,1) = temp(:,:,:,1).*sinTheta;
term2(:,:,:,2) = temp(:,:,:,2).*sinTheta;
term2(:,:,:,3) = temp(:,:,:,3).*sinTheta;
clear temp sinTheta


term3 = zeros(size(v));
temp = sum(axis.*v, 4).*(1-cosTheta);
term3(:,:,:,1) = axis(:,:,:,1).*temp;
term3(:,:,:,2) = axis(:,:,:,2).*temp;
term3(:,:,:,3) = axis(:,:,:,3).*temp;

v_rot = term1 + term2 + term3;

end

