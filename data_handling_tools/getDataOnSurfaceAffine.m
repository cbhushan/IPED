% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [d] = getDataOnSurfaceAffine(m_dfs, m_nii, M_world, origin_loc, s_nii, method)
% Gets data on vertices of m_dfs after transforming it to s_nii volume. 
%

workdir = tempname();
mkdir(workdir);

if ~exist('method', 'var')
   method = 'linear';
end

if ~ischar(m_dfs)
   temp_dfs = fullfile(workdir, [Random_String(16) '.dfs']);
   writedfs(temp_dfs, m_dfs);
   m_dfs = temp_dfs;
end

temp_out_dfs = fullfile(workdir, [Random_String(16) '.dfs']);
dfs_out = affine_transform_dfs(m_dfs, M_world, origin_loc, temp_out_dfs, m_nii, s_nii);
d = getDataOnSurface(dfs_out, s_nii, method);

rmdir(workdir, 's');
end
