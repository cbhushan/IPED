% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function out = copy_nifti_header(header_source, data, type, outputfile)
% Copies parts of nifti header (specified by TYPE) from SOURCE to DATA and
% writes the OUTPUTFILE(optional).
%

workDir = [pwd '/temp_workdir_' Random_String(10)];
mkdir(workDir);

if ischar(header_source)
   header_source = load_untouch_nii_gz(header_source, workDir);
elseif ~isfield(header_source,'untouch') || header_source.untouch ~= 1
   error('Please use ''load_untouch_nii.m'' to load the structure.');
end

if ischar(data)
   data = load_untouch_nii_gz(data, workDir);
elseif ~isfield(data,'untouch') || data.untouch ~= 1
   error('Please use ''load_untouch_nii.m'' to load the structure.');
end

out = data;

if strcmpi(type, 'sform')
   out.hdr.hist.sform_code = header_source.hdr.hist.sform_code;
   out.hdr.hist.srow_x = header_source.hdr.hist.srow_x;
   out.hdr.hist.srow_y = header_source.hdr.hist.srow_y;
   out.hdr.hist.srow_z = header_source.hdr.hist.srow_z;
else
   error(['Unsupported type: ' type])
end

out.hdr.dime.dim = header_source.hdr.dime.dim;
out.hdr.dime.pixdim = header_source.hdr.dime.pixdim;
out.hdr.dime.scl_slope = header_source.hdr.dime.scl_slope;
out.hdr.dime.scl_inter = header_source.hdr.dime.scl_inter;

if exist('outputfile','var')
   save_untouch_nii_gz(out, outputfile, workDir)
end

rmdir(workDir, 's');
end
