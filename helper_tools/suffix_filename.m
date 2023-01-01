% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function out_fname = suffix_filename(fname, suffix)
% Puts suffix at relevant location. useful for filenames with special meaning.

if strcmp(fname(end-11:end), '.mask.nii.gz')
   out_fname = [fname(1:end-12) suffix '.mask.nii.gz'];
   
elseif strcmp(fname(end-12:end), '.label.nii.gz')
   out_fname = [fname(1:end-13) suffix '.label.nii.gz'];
   
elseif strcmp(fname(end-10:end), '.eig.nii.gz')
   out_fname = [fname(1:end-11) suffix '.eig.nii.gz'];

else
   [f, e] = remove_extension(fname);
   out_fname = [f suffix e];
   
end
