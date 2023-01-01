% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function str_out = escape_filename(filename)
% Replaces all occurance of \ by \\ in input string or cellstr.
% It is useful for escaping windows filenames and paths so that fprintf does not generate:
%      Warning: Escape sequence 'U' is not valid
%

if ischar(filename)
   str_out = strrep(filename, '\\', '\'); % to avoid repeating a previously escaped \
   str_out = strrep(str_out, '\', '\\');
   
elseif iscellstr(filename)
   cfuna = @(s) strrep(s, '\\', '\'); % to avoid repeating a previously escaped \
   str_out = cellfun(cfuna, filename, 'UniformOutput', false);
   cfun = @(s) strrep(s, '\', '\\');
   str_out = cellfun(cfun, str_out, 'UniformOutput', false);
else
   error('Input type must be string or cellstring.')
end

end
