% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [nii, reorient_matrix, sform_new, niifileFlag] = load_nii_BIG_Lab(fname)
% Wrapper around many custom functions to ease life. It does following
% while loading nifti-1 file : 
%   1. Adds sform matrix, if sform does not exist.
%   2. Applies reorientation to volume to make sform cannonical (sform & other nifti headers are
%      also appropriately corrected).
%   3. Updates/replaces qform to reflect new reoriented sform.
%
% Details: 
%   fname - filename of nifti file; .nii or .nii.gz
%   nii - loaded nifti structure
%   reorient_matrix, sform_new, niifileFlag - see corresponding functions for details. 
%
% The loaded nifti structure can be saved to disk by calling either of following: 
%   save_untouch_nii_gz(nii, 'output_filename.nii.gz')
%   save_untouch_nii(nii, 'output_filename.nii') % using Jimmy Shen's toolbox
%
% Note that this function: 
%   1. Prefers sform over qform, when both are present (can be modified easily, if required)
%   2. Uses some intermediate file I/O for some operation - so could be a little slower
%

workdir = tempname();
mkdir(workdir);
temp_fname = fullfile(workdir, [randstr(16) '.nii.gz']);

[niifileFlag, outFile] = check_nifti_file(fname, workdir);
[nii, reorient_matrix, sform_new] = reorient_nifti_sform(outFile, temp_fname);

rmdir(workdir, 's');
end
