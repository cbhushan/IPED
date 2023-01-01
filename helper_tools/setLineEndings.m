% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function setLineEndings(filename, lineEnding)
% Sets the line endings of a text file to desired line-ending.
% NOTE: This function OVERWRITES the input file.
%
% Usage:
%   setLineEndings(filename) % uses current OS to infer lineEnding
%   setLineEndings(filename, lineEnding)
%
% Valid values of lineEndings:
%   - 'CRLF' or 'win': Windows line endings (\r\n)
%   - 'LF' or 'unix' : Unix line endings (\n) for Mac and Linux

fid = fopen(filename, 'r');
if (fid == -1)
   fprintf('Could not open file %s. Make sure that it is a valid file.\n', escape_filename(filename));
   return;
end

if nargin<2
   if ispc()
      lineEnding = 'win';
   else % unix
      lineEnding = 'unix';
   end
end

lineEnding = lower(lineEnding);
if ismember(lineEnding, {'crlf', 'win'})
   lineEnding = sprintf('\r\n');
else
   lineEnding = sprintf('\n');
end

fileContents = '';
line = fgetl(fid);
while (ischar(line))
   fileContents = [fileContents line lineEnding];
   line = fgetl(fid);
end
fclose(fid);

% overwrite contents
fid = fopen(filename, 'w');
fprintf(fid, fileContents);
fclose(fid);

end
