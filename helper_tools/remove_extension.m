% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [fileName, ext] = remove_extension( fileName )
% Returns file name (with path) removing following extensions:
%

ext = '';

if length(fileName)>=3 && strcmpi(fileName(end-2:end), '.gz')
   ext = [fileName(end-2:end) ext];
   fileName = fileName(1:end-3);   
end

while length(fileName)>=4 && ( ...
      strcmpi(fileName(end-3:end), '.nii') || strcmpi(fileName(end-3:end), '.img') || ...
      strcmpi(fileName(end-3:end), '.hdr') || strcmpi(fileName(end-3:end), '.dfs') ||...
      strcmpi(fileName(end-3:end), '.dfc') || strcmpi(fileName(end-3:end), '.txt') ||...
      strcmpi(fileName(end-3:end), '.eig') || strcmpi(fileName(end-3:end), '.ext') ||...
      strcmpi(fileName(end-3:end), '.gii') || strcmpi(fileName(end-3:end), '.png'));
   
   ext = [fileName(end-3:end) ext];
   fileName = fileName(1:end-4);
end

if length(fileName)>=5 && strcmpi(fileName(end-4:end), '.bmat')
   ext = [fileName(end-4:end) ext];
   fileName = fileName(1:end-5);
end

end
