% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function flag_cell = readFlagConfigFile(fname)
% Reads BDP flag configuration file and return the flag and parameters in a cell of strings. No
% checks. The file should follow following format: 
% 
%      # Any line starting with # is comment
%      <flag1> [<flag1 parameter>]
%      <flag2> [<flag2 parameter>] # rest of this line is also comment
%      <flagN> [<flagN parameter>]#This is also comment
%

if exist(fname, 'file')~=2   
   error('BDP:FileDoesNotExist', bdp_linewrap(['BDP could not find specified --flag-conf-file. Check to make sure '...
      'that the following file exits: ' escape_filename(fname)]));
end

% dirty workaround different line-endings
flag_str = '';
fid = fopen(fname, 'r');
tline = fgetl(fid);
while ischar(tline)
    flag_str = sprintf('%s\n%s', flag_str, tline);
    tline = fgetl(fid);
end
flag_str = sprintf('%s \n', flag_str);
fclose(fid);

% find comments & throw them away 
% s = regexprep(flag_str, '\s+#[^\n]*\n', '\n'); % shell script style comment rule
s = regexprep(flag_str, '#[^\n]*\n', '\n');

% separate flags
flag_cell = regexpi(s, '\s*(\S+)', 'tokens');

end
