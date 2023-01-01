% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [dwi_out_file, bmat_file, bvec_file, nii_RAS, bMatrice_RAS, bvec_RAS] = reorientDWI_BDP(dwi_file, dwi_out_file, bmat, bdp_opts)
% Same as reorient_nifti_sform followed by bvec rotation, if any. 

fprintf('\nChecking orientation information...')
[nii_RAS, reorient_matrix, sform_new, bMatrice_RAS] = reorient_nifti_sform(dwi_file, dwi_out_file, bmat);

% reorient bvec, if it exists
if isempty(bdp_opts.bvec_file)
   bvec_file = '';
   bvec_RAS = [];
else
   bvec = readBvecBval(bdp_opts.bvec_file);
   bvec_RAS = transpose(reorient_matrix*(bvec'));
   bvec_file = [remove_extension(dwi_out_file) '.bvec'];
   writeBvecBvalFile(bvec_RAS, bvec_file);
end
bmat_file = [remove_extension(dwi_out_file) '.bmat'];

fprintf('Done');
end
