% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function h = displayDiffusionDirections(vec, varargin)
% Sphere plot showing diffusion directions. This function automatically adds diametrically opposite
% points. 
% Usage: 
%    h = displayDiffusionDirections(vec)
%    h = displayDiffusionDirections(vec, markercolor)
%    h = displayDiffusionDirections(vec, markercolor, h) % does not plot sphere. Assumes that a sphere
%                                                     % is already plotted in figure with handle h
%

vec = [vec; -1*vec]; % add diametrically opposite points
h = displayDirections(vec, varargin{:});

end
