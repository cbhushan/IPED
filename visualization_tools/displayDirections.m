% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function h = displayDirections(vec, markercolor, h)
% Plots directions as points on sphere. Also see displayDiffusionDirections.m
% Usage: 
%    h = displayDirections(vec)
%    h = displayDirections(vec, markercolor)
%    h = displayDirections(vec, markercolor, h) % does not plot sphere. Assumes that a sphere
%                                               % is already plotted in figure with handle h
%

marker_size = 7;

if ~exist('markercolor', 'var')
   markercolor = [0 1 0];
end

% draw sphere
if nargin<3
   mat_p = [0.5 0.7 0.25]; % material properties
   [xr, yr, zr] = sphere(30);
   Cm = 0.9*ones(size(xr));
   h = figure;
   surf(xr, yr, zr, Cm,'EdgeColor',[0 0 0], 'LineStyle',':', 'LineWidth',0.1,...
      'FaceLighting','gouraud', 'BackFaceLighting','lit', ...
      'AmbientStrength', mat_p(1), 'DiffuseStrength', mat_p(2), 'SpecularStrength', mat_p(3), ...
      'FaceAlpha', 1);
   colormap gray
   light('Position',[1 -1 3], 'Style','infinite')
end

ax = get(h,'CurrentAxes');
hold(ax, 'on');
plot3(vec(:,1), vec(:,2), vec(:,3), 'LineStyle' , 'none', ...
   'Marker', 'o', 'MarkerSize', marker_size, 'MarkerFaceColor', markercolor, 'MarkerEdgeColor', markercolor)
axis equal
axis vis3d
view([45 45]);
axis off
set(h, 'color', [0 0 0])

end
