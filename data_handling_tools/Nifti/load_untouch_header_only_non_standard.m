% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


% Same as load_untouch_header_only from the nifti toolbox, with line 81
% modified to not throw an error when encounter a non-standard nifti format
% header. Should produce correct header for standard nifti files.
%
function [hdr, ext, filetype, machine] = load_untouch_header_only(filename)

   if ~exist('filename','var')
      error('Usage:  [header, ext, filetype, machine] = load_untouch_header_only(filename)');
   end

   %  Read the dataset header
   %
   [hdr, filetype, fileprefix, machine] = load_nii_hdr(filename);

   if filetype == 0
      hdr = load_untouch0_nii_hdr(fileprefix, machine);
      ext = [];
   else
      hdr = load_untouch_nii_hdr(fileprefix, machine, filetype);

      %  Read the header extension
      %
      ext = load_nii_ext(filename);
   end

   %  Set bitpix according to datatype
   %
   %  /*Acceptable values for datatype are*/ 
   %
   %     0 None                     (Unknown bit per voxel) % DT_NONE, DT_UNKNOWN 
   %     1 Binary                         (ubit1, bitpix=1) % DT_BINARY 
   %     2 Unsigned char         (uchar or uint8, bitpix=8) % DT_UINT8, NIFTI_TYPE_UINT8 
   %     4 Signed short                  (int16, bitpix=16) % DT_INT16, NIFTI_TYPE_INT16 
   %     8 Signed integer                (int32, bitpix=32) % DT_INT32, NIFTI_TYPE_INT32 
   %    16 Floating point    (single or float32, bitpix=32) % DT_FLOAT32, NIFTI_TYPE_FLOAT32 
   %    32 Complex, 2 float32      (Use float32, bitpix=64) % DT_COMPLEX64, NIFTI_TYPE_COMPLEX64
   %    64 Double precision  (double or float64, bitpix=64) % DT_FLOAT64, NIFTI_TYPE_FLOAT64 
   %   128 uint8 RGB                 (Use uint8, bitpix=24) % DT_RGB24, NIFTI_TYPE_RGB24 
   %   256 Signed char            (schar or int8, bitpix=8) % DT_INT8, NIFTI_TYPE_INT8 
   %   511 Single RGB              (Use float32, bitpix=96) % DT_RGB96, NIFTI_TYPE_RGB96
   %   512 Unsigned short               (uint16, bitpix=16) % DT_UNINT16, NIFTI_TYPE_UNINT16 
   %   768 Unsigned integer             (uint32, bitpix=32) % DT_UNINT32, NIFTI_TYPE_UNINT32 
   %  1024 Signed long long              (int64, bitpix=64) % DT_INT64, NIFTI_TYPE_INT64
   %  1280 Unsigned long long           (uint64, bitpix=64) % DT_UINT64, NIFTI_TYPE_UINT64 
   %  1536 Long double, float128  (Unsupported, bitpix=128) % DT_FLOAT128, NIFTI_TYPE_FLOAT128 
   %  1792 Complex128, 2 float64  (Use float64, bitpix=128) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
   %  2048 Complex256, 2 float128 (Unsupported, bitpix=256) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
   %
   switch hdr.dime.datatype
   case   1,
      hdr.dime.bitpix = 1;  precision = 'ubit1';
   case   2,
      hdr.dime.bitpix = 8;  precision = 'uint8';
   case   4,
      hdr.dime.bitpix = 16; precision = 'int16';
   case   8,
      hdr.dime.bitpix = 32; precision = 'int32';
   case  16,
      hdr.dime.bitpix = 32; precision = 'float32';
   case  32,
      hdr.dime.bitpix = 64; precision = 'float32';
   case  64,
      hdr.dime.bitpix = 64; precision = 'float64';
   case 128,
      hdr.dime.bitpix = 24; precision = 'uint8';
   case 256 
      hdr.dime.bitpix = 8;  precision = 'int8';
   case 511 
      hdr.dime.bitpix = 96; precision = 'float32';
   case 512 
      hdr.dime.bitpix = 16; precision = 'uint16';
   case 768 
      hdr.dime.bitpix = 32; precision = 'uint32';
   case 1024
      hdr.dime.bitpix = 64; precision = 'int64';
   case 1280
      hdr.dime.bitpix = 64; precision = 'uint64';
   case 1792,
      hdr.dime.bitpix = 128; precision = 'float64';
   otherwise
      warning('This is NOT a standard nifti file.\n\n    hdr.dime.datatype = %d \n\nLoaded header might have some problems. \n',hdr.dime.datatype); 
   end

   tmp = hdr.dime.dim(2:end);
   tmp(find(tmp < 1)) = 1;
   hdr.dime.dim(2:end) = tmp;


   return					% load_untouch_header_only
