% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [pathstr, baseName] = fileBaseName( fileName )
% Returns path string & the base name of the fileName (removing multiple
% extensions - see remove_extension.m)
%
% Usage:
%     [pathstr, baseName ] = file_base_name(fileName)
%     baseName = file_base_name(fileName)
%


[pathstr, baseName, ext] = fileparts(fileName);
baseName = remove_extension([baseName ext]);

if (nargout <= 1)
   pathstr = baseName;
end


end

