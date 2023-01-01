% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function img_out = imgaussian3D(img, sgma, ksiz)
% A simple wrapper around imfilter for 3D Gasussian smoothing. Uses separable property of
% Gaussian for faster computation. 
%  
% This function is slower than imgaussian() mex code, which produced identical results (under
% numerical errors of 1e-12). 
%

if length(sgma)~=1
   error('sgma must be scalar - Gaussian is (trivially) separable only for isotropic smoothing.')
end

if ~exist('siz','var')
   ksiz = sgma*6;
end

x = -ceil(ksiz/2):ceil(ksiz/2);
H = exp(-(x.^2/(2*sgma^2)));
H = H/sum(H(:));

Hx = reshape(H,[length(H) 1 1]);
Hy = reshape(H,[1 length(H) 1]);
Hz = reshape(H,[1 1 length(H)]);

img_out = imfilter(imfilter(imfilter(img, Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate'),Hz, 'same' ,'replicate');

end
