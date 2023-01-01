% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function v_conn = neighborsVerticesCommonFaces(dfs, N)
% Returns neighboring vertices corresponding to each vertex in dfs. Two
% vertices are considered as neighbors if they share an edge. 
% N defines order/length of connectivity - N=1 one is same as immediate
%   neighbors; N=2 includes neighbors of neighbors; and likewise
%
% Open matlabpool for speed.
%

if ischar(dfs)
   dfs = readdfsGz(dfs);
elseif ~isfield(dfs,'vertices') || ~isfield(dfs,'faces')
   error('dfs structure must have vertices and faces field!');
end

n_vert = length(dfs.vertices);
r = [dfs.faces(:,1);dfs.faces(:,1); dfs.faces(:,2);dfs.faces(:,2); dfs.faces(:,3);dfs.faces(:,3)];
c = [dfs.faces(:,2);dfs.faces(:,3); dfs.faces(:,1);dfs.faces(:,3); dfs.faces(:,1);dfs.faces(:,2)];
d = ones(size(r));
C = sparse(r, c, d, n_vert, n_vert); % symmetric matrix
C = spones(C);

[r, c, d] = find(C);
if length(unique(c))~=n_vert
   error('Isolated vertices found - this implementation wont work! May be first run my_clean_patch?')
end


b_ind = [0; find(diff([c; 0]))];
v_conn = cell(n_vert,1); % one empty cell per vertex
for i = 1:n_vert
   v_conn{i} = [i; r((b_ind(i)+1) : b_ind(i+1))]; % add itself to connectivity matrix
end

% Order N connectivity
v_conn1 = v_conn; % 1st neighbors
for n = 2:N
   v_conn_update = cell(n_vert,1); % one empty cell per vertex
   parfor i = 1:n_vert
      v_conn_update{i} = unique(vertcat(v_conn1{v_conn{i}}));
   end
   v_conn = v_conn_update;
end

end

