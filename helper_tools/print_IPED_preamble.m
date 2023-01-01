% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function print_IPED_preamble(outFile)
% Print short intro with license

if exist('outFile', 'var')
    try
        fid = fopen(outFile, 'w');
    catch err
        error('BDP:FileDoesNotExist', 'Could not create %s for writing preamble', outFile);
    end
else
    fid = 1; % No file given, so print to screen
end

% show short license content
[h_dir, ff, ee] = fileparts(mfilename('fullpath'));
short_lic_file = fullfile(h_dir, 'license_short.txt');

wdt = 80;
fin = fopen(short_lic_file, 'r');
lic_s = '';
while ~feof(fin)
   tline = fgetl(fin);
   wdt = max(wdt, length(tline));
   lic_s = [lic_s tline '\n'];
end
fclose(fin);

sep_str(1:wdt) = '=';
fprintf(fid, '\n%s\n', sep_str);
fprintf(fid, lic_s);
fprintf(fid, '\n%s\n', sep_str);

