% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function bdpPrintWarning(hdr, msg, wdt)
% Prints warning on command line outputs.

if ~exist('wdt', 'var')
   wdt = 80;
end
sep_str(1:wdt) = '*';

hdr = ['WARNING : ' hdr];
sz = length(hdr);

if wdt>sz
   blk(1:ceil((wdt-sz)/2)) = ' ';
else
   blk = '';
end
hdr = ['*' blk hdr '\n'];
strout = bdp_linewrap(['*  ' msg], wdt-5, '\n*  ');

fprintf('\n%s\n', sep_str);
fprintf(hdr);
fprintf('*\n');
fprintf(strout);
fprintf('\n%s\n', sep_str);

end
