% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [T, d_out] = vertices2faces(dfs, data)
% Returns a matrix operator which transfers data from vertex of triangulated mesh to faces. Data on
% a face is average of the data on its vertices.

if ischar(dfs)
   dfs = readdfsGz(dfs);
elseif ~isfield(dfs,'vertices') || ~isfield(dfs,'faces')
   error('dfs structure must have vertices and faces field!');
end

if exist('data','var') && (length(dfs.vertices)~=length(data))
   error('data must be of same length as number of vertices.');
end

nF = length(dfs.faces);
i = [1:nF 1:nF 1:nF]';
j = dfs.faces(:);
s = ones(size(i))/3;
T = sparse(i, j, s, nF, length(dfs.vertices));

if exist('data','var')
   d_out = T*data;
end

end
