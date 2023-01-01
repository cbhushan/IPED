% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [T, d_out] = faces2vertices(dfs, data)
% Returns a matrix operator which transfers data from faces of triangulated mesh to vertices. Data
% on a vertex is average of all faces attached to the vertex. 

if ischar(dfs)
   dfs = readdfsGz(dfs);
elseif ~isfield(dfs,'vertices') || ~isfield(dfs,'faces')
   error('dfs structure must have vertices and faces field!');
end

if exist('data','var') && (length(dfs.faces)~=length(data))
   error('data must be of same length as number of vertices.');
end

nF = length(dfs.faces);
nV = length(dfs.vertices);

TR = triangulation(dfs.faces, dfs.vertices);
clear dfs
ti = vertexAttachments(TR);

[~, nc] = cellfun(@size, ti);
i = cell(nV, 1);
s = cell(nV, 1);
for v = 1:nV
   i{v} = ones(1,nc(v))*v;
   s{v} = ones(1,nc(v))/nc(v);
end

T = sparse([i{:}], [ti{:}], [s{:}], nV, nF);

if exist('data','var')
   d_out = T*data;
end

end
