% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function dataMasked = mask_nii(dataIn, maskIn, fileName, gz)
%MASK_NII Masks the nii image using a binary/logical mask. Accepts both nii
%structure and .nii or .nii.gz files. writes the output to disk.
%
%   dataIn - nii structure or filename (.nii or .nii.gz) with path for nii
%            image to be masked 
%
%   maskIn - nii structure or filename (.nii or .nii.gz) with path for the
%            mask (binary)
%
%   gz - set to false to disable compressesd output
%
%   dataMasked - nii structure of the masked image
%

workDir = tempname; %[pwd '/temp_workdir_' Random_String(10)];
mkdir(workDir);


if ischar(dataIn)
   data = load_untouch_nii_gz(dataIn, workDir);
else
   if ~isfield(dataIn,'untouch') || dataIn.untouch ~= 1
      error('Please use ''load_untouch_nii.m'' to load the file.');
   end
   data = dataIn;
   clear dataIn
end


if ischar(maskIn)
   mask = load_untouch_nii_gz(maskIn, workDir);
else
   if ~isfield(maskIn,'untouch') || maskIn.untouch ~= 1
      error('Please use ''load_untouch_nii.m'' to load the file.');
   end
   mask = maskIn;
   clear maskIn
end


dataMasked = data;

dataDim = size(data.img);
maskDim = size(mask.img);

dataMasked.img = data.img;
if length(dataDim)==length(maskDim) && isequal(dataDim, maskDim)
  dataMasked.img(mask.img<=0) = 0;
  
elseif ((length(dataDim)-length(maskDim))==1 && isequal(dataDim(1:end-1), maskDim))
  maskImg = mask.img>0;
  maskImg = repmat(maskImg, [ones([1 length(maskDim)]) dataDim(end)]);
  dataMasked.img(~maskImg)=0;
  
else
  error('Dimension mismatch.');
end

if exist('fileName','var')
    if exist('gz', 'var') && ~gz
      save_untouch_nii_wrapper(dataMasked, fileName);
   else
      save_untouch_nii_gz(dataMasked, fileName, workDir);
    end
end
  
rmdir(workDir, 's');

end

