% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [ savedFileName] = save_nii_gz( niiData, fileName, datatype )
% Wrapper around save_nii to create compressed GNU zip files. Only works
% with nii structures.
%
%   niiData: nii structure or .nii filename with full path.
%
%   fileName: Destination filename, with/without '.nii' or '.nii.gz'.
%             Use [] for no filename. When [] is used, filename is
%             retrived from the nii structure or the file is compressed
%             with the same filename (as that of .nii).
%
%   datatype: [Optional] NiFTI datatype; must be integer. (see save_nii.m
%             for exhaustive list) Few common data types:
%                4    Signed short; int16; NIFTI_TYPE_INT16 (usually for masks)
%               16    Floating point; single or float32; NIFTI_TYPE_FLOAT32
%               64    Double precision; double or float64; DT_FLOAT64, NIFTI_TYPE_FLOAT64
%

if isstruct(niiData)
   
   if isempty(fileName)
      fileName = niiData.fileprefix;
   else
      fileName = remove_extension(fileName);
      niiData.fileprefix = fileName;
   end
   
   niiData.hdr.dime.cal_max = max(niiData.img(:));
   niiData.hdr.dime.cal_min = min(niiData.img(:));
   
   if exist('datatype','var')
      if (niiData.hdr.dime.datatype == 128) && (datatype~=128)
         niiData.hdr.dime.dim(5) = 3;
         niiData.hdr.dime.dim(1) = 4;
      end
      niiData.hdr.dime.datatype = datatype;
      niiData.hdr.dime.bitpix = datatype2bitpix(datatype);
   end
   
   save_nii_wrapper(niiData, [fileName '.nii']);
   fname = gzip([fileName '.nii']);
   savedFileName = fname{1};
   delete([fileName '.nii'])
   
else
   error('niiData must be nifti structure')
end

end

