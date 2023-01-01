% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function D = createDonSurfaceVertexApprox(surf)
% Creates an approx. derivative operator on surface vertices (for data defined on vertices). This
% just a proxy and is NOT accurate derivative. Could be useful for finding edges of surface labels.
% Computes following for each vertex:
%
%      Di = Sum over neighbors (Vi - Vn)
%
%    where, Vi = function value at vertex i
%           Vn = function value at a neighbor of vertex i
%           Di = sum of pair-wise difference with neighbors at each vertex

n_vert = size(surf.vertices, 1);

r = [surf.faces(:,1);surf.faces(:,1); surf.faces(:,2);surf.faces(:,2); surf.faces(:,3);surf.faces(:,3)];
c = [surf.faces(:,2);surf.faces(:,3); surf.faces(:,1);surf.faces(:,3); surf.faces(:,1);surf.faces(:,2)];
d = ones(size(r));
C = sparse(r, c, d, n_vert, n_vert); % symmetric matrix
C = spones(C); % connectivity at each vertex (except itself)

n_conn = full(sum(C,2));

D = spdiags(n_conn, 0, n_vert, n_vert) - C;

end
