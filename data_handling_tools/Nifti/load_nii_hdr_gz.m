% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [ data ] = load_nii_hdr_gz( file )
%wrapper around load_nii_hdr to work with .gz files.
%
%   1. Extracts GNU zip files, if file is zipped.
%   2. Loads the data & returns it.
%   3. Deletes the unzipped file, if file is zipped.
%
%   file - is a string with relative / absolute path to .gz or .nii file
%
%   Requires NIFTI toolbox.

%fprintf('\nLoading file:\n %s\n',file);

loc = strfind(file, '.gz');
len = length(file);
if (~isempty(loc)) && (len-loc(end) == 2)  % .gz extension found
  niiData = gunzip(file,'temp');
  data = load_nii_hdr(char(niiData));
  delete(char(niiData));
  data.fileprefix = data.fileprefix(6:end); % remove 'temp/'
else
  data = load_nii_hdr(file);
end

%fprintf('Done...\n');
end

