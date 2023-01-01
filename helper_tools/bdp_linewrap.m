% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function strout = bdp_linewrap(str, maxchars, separator)
% Wrapper around linewrap and strjoin to simplify linewrap

if nargin == 2 && ischar(maxchars)
   separator = maxchars;
   clear maxchars
end

if ~exist('maxchars', 'var')
   maxchars = 80;
end

if ~exist('separator', 'var')
   separator = '\n';
end

if ischar(str)
   strout = strjoin_KY(linewrap(str, maxchars), separator);
   
elseif iscellstr(str)
   strout = cellfun(@(x) bdp_linewrap(x, maxchars, separator), str, 'UniformOutput', false);
   strout = strjoin_KY(strout, '');
   
else
   error('Input str must be string or cellstr.')
end

end
