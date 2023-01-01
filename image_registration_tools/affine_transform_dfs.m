% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function dfs_out = affine_transform_dfs(moving_dfs_file, F, origin_loc, output_dfs_file, moving_nii, static_nii, gz)
% This function (affine) transforms the MOVING_DFS_FILE (corresponding to MOVING_NII) to
% coordinate space of STATIC_NII such that output surface overlays correctly in BrainSuite.
%    F - 4x4 affine matrix transformation
%    origin_loc - 3x1 vector representing origin. 
%
% See affine_transform_nii.m for important note about [F, origin]. It would be typically output of
% rigid registration between MOVING_NII and STATIC_NII (register_file_rigid.m). Following equations
% describe the relation between BrainSuite's coordinate system and NifTI coordinates. 
%
%         Xw = Sn * Xn
%            = Mrest*Sc * Xn
%            = Mrest*(Sc*Xn)
%            = Mrest * Xb
%
%     Xw - world coordinate (nifti coordinate) of moving_nii
%     Xn - voxel coordinate in RAS configurate (as BrainSuite operates in RAS) of moving_nii
%          Indexing starts from zero (0)
%     Sn - sform transformation in nifti header corresponding to RAS configuration of voxels of
%          moving_nii (see reorient_nifti_sform.m to obtain Sn, Xn for an arbitrary nifti file)
%     Sc - Scaling component of Sn
%     Mrest - Residual transformation after taking out scaling component(Sc)
%     Xb - BrainSuite's coordinate (for surfaces) of moving_nii
%
% Now, after applying affine transformation the coordinates should match to that of static file i.e.
% Xw_m = Xw_s. Also see affine_transform_nii.m 
%
%                Xw_s = Xw_m
%                Xw_s = (Sinv*F*S) * Xw
%                Xw_s = (Sinv*F*S) * Mrest * Xb
%      Mrest_s * Xb_s = (Sinv*F*S) * Mrest * Xb
%                Xb_s = inv(Mrest_s) * (Sinv*F*S) * Mrest * Xb
%                Xb_s = Mfinal * Xb
%
%     Xw_m - world coordinate of moving_nii after affine transformation 
%     Xw_s - world coordinate of static_nii; Xw_s = Mrest_s * Xb_s
%     S  - translation matrix (4x4) which transfers origin to ORIGIN_LOC
%     F  - Affine transform to be applied
%     Sinv - inv(S); Puts back the origin from ORIGIN_LOC to [0; 0; 0]
%     Xb_s - BrainSuite's coordinate (for surfaces) of static_nii
%     Mrest_s - Equivalent of Mrest but for static_nii
%

workdir = tempname();
mkdir(workdir);
tempfile = fullfile(workdir, [Random_String(16) '.nii.gz']);

if ~exist('gz', 'var')
   gz = false;
end

if ~ischar(moving_nii)
   moving_nii = save_untouch_nii_wrapper(moving_nii, fullfile(workdir, [Random_String(16) '.nii.gz']));
end

if ~ischar(static_nii)
   static_nii = save_untouch_nii_wrapper(static_nii, fullfile(workdir, [Random_String(16) '.nii.gz']));
end

[nii_mv, ~, ST_mv] = reorient_nifti_sform(moving_nii, tempfile);
[nii_st, ~, ST_st] = reorient_nifti_sform(static_nii, tempfile);

Mrest_mv = ST_mv*diag([1./abs(nii_mv.hdr.dime.pixdim(2:4)) 1]);
Mrest_st = ST_st*diag([1./abs(nii_st.hdr.dime.pixdim(2:4)) 1]);

S = eye(4);
S(1:3, 4) = -1* origin_loc;

dfs_in = readdfsGz(moving_dfs_file);
vert_mv = dfs_in.vertices;
vert_mv(:,4) = 1;
vert_mv = vert_mv';

M_final = inv(Mrest_st)*inv(S)*F*S*Mrest_mv;
vert_st = M_final*vert_mv;

dfs_out = dfs_in;
dfs_out.vertices = [vert_st(1:3,:)]';

if isfield(dfs_in, 'normals')
   normals = dfs_in.normals;
   normals(:,4) = 1;
   normals = normals';
   normals_transform = M_final*normals;
   dfs_out.normals = [normals_transform(1:3,:)]';
end

if gz % gzipped output
   writedfsGz(output_dfs_file, dfs_out);
else
   writedfs(output_dfs_file, dfs_out);
end

rmdir(workdir, 's');
clearvars -except dfs_out

end
