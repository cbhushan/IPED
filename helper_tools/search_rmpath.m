% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function dir_fullpath = search_rmpath(dir_name)
% Searches for dir_name in matlab path and if found, removes dir_name and
% all its subdirectories from the MATLAB search path.
% Returns the list of directories which were removed. The directories can
% be added to path again by running:
%        addpath(dir_fullpath)
%


% get path as long string
p=path;

f = strfind(p, [dir_name pathsep]);
dir_fullpath = [];

if isempty(f)
   return
   %fprintf('%s - not found in path.\n', dir_name);
else
   % divide string to directory names
   delim=[0 strfind(p, pathsep) length(p)+1];

   for nf = 1:length(f)
      n = 1;
      while delim(n)<f(nf)
         n = n+1;
      end
      dir_root_full = p(delim(n-1)+1:delim(n)-1); % full path to dir_name

      for i=max(n-1, 2):length(delim)
         direc = p(delim(i-1)+1:delim(i)-1);
         if strncmpi(direc, dir_root_full, length(dir_root_full))   
            dir_fullpath = [dir_fullpath direc pathsep];
            rmpath(direc);
         end
      end
   end
   dir_fullpath(end) = []; % remove last pathsep
end

end
