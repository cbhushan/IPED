% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [ data ] = load_nii_gz(file, workdir)
% wrapper around load_nii_gz to work with .gz files.
%   1. Extracts GNU zip files, if file is zipped.
%   2. Loads the data & returns it.
%   3. Deletes the unzipped file, if file is zipped.
%
%   file - is a string with relative / absolute path to .gz or .nii file
%
%   Requires NIFTI toolbox.


loc = strfind(file, '.gz');
len = length(file);
if (~isempty(loc)) && (len-loc(end) == 2)  % .gz extension found
   
   if ~exist('workdir','var')
      workdir = [pwd '/temp_workdir_' Random_String(10)];
      temp_folder = 1;
   end
   
   niiFile = gunzip(file, workdir);
   data = load_nii(char(niiFile));
   delete(char(niiFile));
   
   % fix data.fileprefix to reflect the actual folder
   data.fileprefix = remove_extension(file);
   
   if exist('temp_folder','var') && temp_folder == 1 % delete temp_dir if created
      rmdir(workdir, 's')
   end
   
else
   data = load_nii(file);
end


end

