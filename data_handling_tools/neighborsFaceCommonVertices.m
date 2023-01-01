% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function N = neighborsFaceCommonVertices(dfs)
% Returns neighboring faces corresponding to each face in dfs. Two faces are considered as
% neighbors if they share one or more vertices. 
%

if ischar(dfs)
   dfs = readdfsGz(dfs);
elseif ~isfield(dfs,'vertices') || ~isfield(dfs,'faces')
   error('dfs structure must have vertices and faces field!');
end

TR = triangulation(dfs.faces, dfs.vertices);
ti = vertexAttachments(TR);
tif = [ti(dfs.faces(:,1)) ti(dfs.faces(:,2)) ti(dfs.faces(:,3))];

% Is this part of code vectorizable? 
N = cell(size(dfs.faces,1),1);
for iface = 1:size(dfs.faces,1)
   temp = unique(horzcat(tif{iface,:}));
   N{iface} = temp(temp~=iface); % remove the current face itself
end

end

