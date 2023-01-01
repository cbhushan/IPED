% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [Dx_local, Dy_local, x_unit, y_unit] = createDonSurfaceFace(surf)
% Creates derivative operators, which computes derivative of functions defined at each vertices. The
% derivate is computed w.r.t. a local coordinate system defined at each face of the surface mesh and
% the local coordinate system is different for each face. These derivatives should be first
% expressed in world coordinates before any one-to-one comparison, as shown below. 
%
%   Du_x_local = Dx_local * surf.u(:); % local x-derivative of surf.u
%   Du_y_local = Dy_local * surf.u(:); % local y-derivative of surf.u
%
%   % Derivative of surf.u expressed in world coordinates
%   Du_world = (x_unit.* Du_x_local(:, [1 1 1])) + (y_unit .* Du_y_local(:, [1 1 1])); 
%
% Assumptions: The function is assumed to piece-wise linear over each traingle/face. Hence, the
%              derivative for each face on the surface mesh is constant.
%

X = surf.vertices(:,1);
Y = surf.vertices(:,2);
Z = surf.vertices(:,3);
NumFaces = size(surf.faces,1);
NumVertices = size(X,1);
vertx_1 = surf.faces(:,1); % vertex-index of first vertex of face
vertx_2 = surf.faces(:,2); % vertex-index of second vertex of face
vertx_3 = surf.faces(:,3);
V1 = [X(vertx_1),Y(vertx_1),Z(vertx_1)]; % Vertices of the first triangle of each face
V2 = [X(vertx_2),Y(vertx_2),Z(vertx_2)]; % Vertices of the 2nd triangle of each face
V3 = [X(vertx_3),Y(vertx_3),Z(vertx_3)]; % Vertices of the 3rd triangle of each face


% Local (2D) coordinate for each triangle/face: x-axis along V1-V2; y-axis along normal to the x-axis
x1 = zeros(NumFaces,1); % local x-coordinate of V1; V1 reprents origin
y1 = zeros(NumFaces,1); % local y-coordinate of V1; V1 reprents origin

v2_v1temp = V2-V1;
x2 = sqrt(sum(v2_v1temp.^2, 2)); % local x-coordinate of V2; same as distance b/w V1 and V2
y2 = zeros(NumFaces,1); % local y-coordinate of V2; local x-axis lies along V1-V2, hence zero y-coordinate
x_unit = v2_v1temp ./ x2(:,[1 1 1]);

v3_v1temp = V3-V1;
x3 = dot(v3_v1temp, x_unit, 2); % local x-coordinate of V3; same as projection of V3 on x-axis
mynorm = cross(x_unit, v3_v1temp, 2);
ydir = cross(mynorm, x_unit, 2);
ymag = sqrt(sum(ydir.^2, 2));
y_unit = ydir ./ ymag(:,[1 1 1]);
y3 = dot(y_unit, v3_v1temp, 2); % local y-coordinate of V3;

% Row 1 - for derivative along local x-coordinate, Joshi et al. 2004 - p-harmonic
DT = abs(x1.*y2 - y1.*x2 + x2.*y3 - y2.*x3 + x3.*y1 - y3.*x1); % 2 * area of triangle
tmp_A = [y2-y3, y3-y1, y1-y2]./DT(:,[1 1 1]);
tmp_A = tmp_A(:);

% Row 2 - for derivative along local y-coordinate
tmp_B = [x3-x2, x1-x3, x2-x1]./DT(:,[1 1 1]);
tmp_B = tmp_B(:);

rowno = 1:NumFaces;
rowno_all = [rowno(:); rowno(:); rowno(:)];
vertx_all=[vertx_1; vertx_2; vertx_3];
Dx_local = sparse(rowno_all, vertx_all, tmp_A, NumFaces, NumVertices); % operator for derivative along local x-coordinate
Dy_local = sparse(rowno_all, vertx_all, tmp_B, NumFaces, NumVertices); % operator for derivative along local y-coordinate

clearvars -except Dx_local Dy_local y_unit x_unit
end

