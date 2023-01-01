% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function bitpix = datatype2bitpix (datatype)
% Returns bitpix for the input datatype according to following table: 
% 
%     0 None                     (Unknown bit per voxel) % DT_NONE, DT_UNKNOWN 
%     1 Binary                         (ubit1, bitpix=1) % DT_BINARY 
%     2 Unsigned char         (uchar or uint8, bitpix=8) % DT_UINT8, NIFTI_TYPE_UINT8 
%     4 Signed short                  (int16, bitpix=16) % DT_INT16, NIFTI_TYPE_INT16 
%     8 Signed integer                (int32, bitpix=32) % DT_INT32, NIFTI_TYPE_INT32 
%    16 Floating point    (single or float32, bitpix=32) % DT_FLOAT32, NIFTI_TYPE_FLOAT32 
%    32 Complex, 2 float32      (Use float32, bitpix=64) % DT_COMPLEX64, NIFTI_TYPE_COMPLEX64
%    64 Double precision  (double or float64, bitpix=64) % DT_FLOAT64, NIFTI_TYPE_FLOAT64 
%   128 Red-Green-Blue            (Use uint8, bitpix=24) % DT_RGB24, NIFTI_TYPE_RGB24 
%   256 Signed char            (schar or int8, bitpix=8) % DT_INT8, NIFTI_TYPE_INT8 
%   512 Unsigned short               (uint16, bitpix=16) % DT_UNINT16, NIFTI_TYPE_UNINT16 
%   768 Unsigned integer             (uint32, bitpix=32) % DT_UNINT32, NIFTI_TYPE_UNINT32 
%  1024 Signed long long              (int64, bitpix=64) % DT_INT64, NIFTI_TYPE_INT64
%  1280 Unsigned long long           (uint64, bitpix=64) % DT_UINT64, NIFTI_TYPE_UINT64 
%  1536 Long double, float128  (Unsupported, bitpix=128) % DT_FLOAT128, NIFTI_TYPE_FLOAT128 
%  1792 Complex128, 2 float64  (Use float64, bitpix=128) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
%  2048 Complex256, 2 float128 (Unsupported, bitpix=256) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
%

if datatype == 0
   bitpix = 0;
elseif datatype == 1
   bitpix =  1;
elseif datatype == 2
   bitpix = 8;
elseif datatype == 4
   bitpix = 16;
elseif datatype == 8
   bitpix = 32;
elseif datatype == 16
   bitpix = 32;
elseif datatype == 32
   bitpix = 64;
elseif datatype == 64
   bitpix = 64;
elseif datatype == 128
   bitpix = 24;
elseif datatype == 256
   bitpix = 8;
elseif datatype == 512
   bitpix = 16;
elseif datatype == 768
   bitpix = 32;
elseif datatype == 1024
   bitpix = 64;
elseif datatype == 1280
   bitpix = 64;
elseif datatype == 1536
   bitpix = 128;
elseif datatype == 1792
   bitpix = 128;
elseif datatype == 2048
   bitpix = 256;
else 
   error('\nUnsupported or unknown datatype. Check save_nii.m for details.\n');
end
   
end
   
   
   
   
