% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [ I, pngFileName ] = vol2png_selective(vol, CLim, slices1, slices2, slices3, pngFileName, overlay)
% Generates png format of the slices.
%   vol - nii structure
%   slices1 - vector containing slices numbers to display along first dimension
%   slices2 - along 2nd dimension
%   slices3 - along 3rd dimension
%   CLim - Vector of length 2 for lower limit & upper limit of
%          intnsity. When not specified uses max and min intensity values
%          as limits.
%   pngFileName - [Optional] Saving location of png with full path. When
%                 not specified saves with name [vol.fileprefix '.png']
%   overlay - [optional] Image (m x n) matrix to overlay on bottom
%             right corner of image. Can be used for adding subject name
%             etc.

fprintf('\n==== Writing png files====\n');

if nargin == 6 && ~ischar(pngFileName)
  overlay = pngFileName;
  clear pngFileName
end

if ~isvector(CLim)
  CLim(1) = min(vol.img(:));
  CLim(2) = max(vol.img(:));
end

volThresh = (vol.img - CLim(1)) / (CLim(2)-CLim(1));

volThresh(volThresh>1) = 1;
volThresh(volThresh<0) = 0;

E = permute(volThresh, [2 1 4 3]);
E = flipdim(E, 1);
I1 = montage_image(E , 'size', [1 NaN], 'Indices', slices3);

E = permute(volThresh, [3 1 4 2]);
E = flipdim(E, 1);
I2 = montage_image(E , 'size', [1 NaN], 'Indices', slices2);

E = permute(volThresh, [3 2 4 1]);
E = flipdim(E, 1);
E = flipdim(E, 2);
I3 = montage_image(E , 'size', [1 NaN], 'Indices', slices1);

m = size(I1,1) + size(I2,1) + size(I3,1) + 4;
n = max([size(I1,2) size(I2,2) size(I3,2)]);

I = zeros([m n]);
I(1:size(I1,1), 1:size(I1,2)) = I1;
I((size(I1,1)+3):(size(I1,1)+2+size(I2,1)), 1:size(I2,2)) = I2;
I(m-size(I3,1)+1:end, 1:size(I3,2)) = I3;

if ~exist('pngFileName', 'var')
  pngFileName = [remove_extension(vol.fileprefix) '.color.png'];
end

if exist('overlay', 'var')
  [p, q, ~] = size(overlay);
  I(m-p+1:m, n-q+1:n) = overlay;
end

imwrite(I, pngFileName, 'png');

fprintf('Done...\n');
end

