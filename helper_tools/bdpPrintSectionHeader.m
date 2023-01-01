% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function bdpPrintSectionHeader(hdr, wdt)
% Prints section header. For BDP command line outputs.

if ~exist('wdt', 'var')
   wdt = 80;
end
sep_str(1:wdt) = '=';
sz = length(hdr);

if sz <= wdt
   blk(1:floor((wdt-sz)/2)) = ' ';
   hdr = [blk hdr '\n'];
   
else % sz>wdt
   k = strfind(hdr(floor(sz/2):end), ' ');
   md = floor(sz/2)+k(1)-2;
   hdr1 = hdr(1:md);
   hdr2 = hdr(md+1:end);   
   
   blk1(1:floor((wdt-length(hdr1))/2)) = ' ';
   blk2(1:floor((wdt-length(hdr2))/2)) = ' ';
   
   hdr = [blk1 hdr1 '\n ' blk2 hdr2 '\n'];   
end

fprintf('\n%s\n', sep_str);
fprintf(hdr);
fprintf('%s\n', sep_str);

end
