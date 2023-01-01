% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function map_final = map_of_map(map1, map2)
% Computes MAP2(MAP1())
%
%  MAP1: A -> S  (S warped to A)
%  MAP2: S -> D  (D warped to S), then
%  MAP_FINAL: A -> D
%
% Assuming that S is standard unit spaced ndgrid type grid point starting
% from 1. All the maps are voxel-to-voxel coordinate map. 
%


xmap = double(squeeze(map1(:,:,:,1))); 
ymap = double(squeeze(map1(:,:,:,2))); 
zmap = double(squeeze(map1(:,:,:,3))); 

% map_final(:,:,:,1) = trilinear(double(squeeze(map2(:,:,:,1))), ymap, xmap, zmap);
% map_final(:,:,:,2) = trilinear(double(squeeze(map2(:,:,:,2))), ymap, xmap, zmap);
% map_final(:,:,:,3) = trilinear(double(squeeze(map2(:,:,:,3))), ymap, xmap, zmap);

map_final(:,:,:,1) = interp3(double(squeeze(map2(:,:,:,1))), ymap, xmap, zmap, 'linear', 0);
map_final(:,:,:,2) = interp3(double(squeeze(map2(:,:,:,2))), ymap, xmap, zmap, 'linear', 0);
map_final(:,:,:,3) = interp3(double(squeeze(map2(:,:,:,3))), ymap, xmap, zmap, 'linear', 0);

map_final(isnan(map_final)) = 0;

end
