% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function estimate_mADC(dwi_file, bMatrices)
% mADC

% File checks
if exist(dwi_file, 'file')~=2
   error('DWI file does not exists: %s', dwi_file);
end

if ischar(bMatrices)
   bMatrices = readBmat(bMatrices);
end
DEout = checkDiffusionEncodingScheme(bMatrices);
dwi = load_untouch_nii_gz(dwi_file, true);

sz = size(dwi.img);
nDir = size(dwi.img, 4);
temp = permute(reshape(double(dwi.img), [], nDir), [2,1]); % 1st dim is different diffusion weighting

adc = dwi; clear dwi
adc.hdr.dime.dim(1) = 3;
adc.hdr.dime.dim(5) = 1;
adc.img = compute_mADC(temp, DEout);
adc.img = reshape(adc.img(:), sz(1:3));

fname = save_untouch_nii_gz(adc, [remove_extension(dwi_file) '.ADC.nii.gz'], 16);
fprintf('\nWritten ADC file: %s\n', fname);
