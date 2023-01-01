% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [ outFile ] = clear_qform( dataFile, outFile )
% Clear qform values (sets all of them to 0) if sform is set in a standard NIFTI
% file. Saves the data in .gz file.
%
%   dataFile - .nii or .nii.gz filename with full path
%   outFile -  output file name.
%

data = load_untouch_nii_gz(dataFile);

if (data.hdr.hist.sform_code <= 0)
   fprintf('\n sform is NOT set. Clearing qform values can make header useless. Terminating without any changes.\n');
else

	data.hdr.hist.qform_code = 0;
	data.hdr.hist.qoffset_x = 0;
	data.hdr.hist.qoffset_y = 0;
	data.hdr.hist.qoffset_z = 0;
	
	data.hdr.hist.quatern_b = 0;
	data.hdr.hist.quatern_c = 0;
	data.hdr.hist.quatern_d = 0;
	
	data.hdr.hist.qform_code = 0;
	data.hdr.dime.pixdim(1) = 0;
	
	save_untouch_nii_gz(data, outFile);
end

end

