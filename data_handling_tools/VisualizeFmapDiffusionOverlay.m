% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [FmapMag_Dcoord, dwi] = VisualizeFmapDiffusionOverlay(dwi_file, bvec, bval, FmapMag_file, FmapPh_file, echo_space, PED, output_filebase)
%
% Saves overlay of distorted fieldmap magnitude image and diffusion images.
%
% echo_space - Echo spacing in sec.
% PED can be one of following: x, x-, y, y-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

workdir = tempname();
mkdir(workdir);
output_filebase = remove_extension(output_filebase);

% Reorientation is imp to interpret PED correctly
FmapMag = reorient_nifti_sform(FmapMag_file, fullfile(workdir, fileBaseName(FmapMag_file)));
FmapPh = reorient_nifti_sform(FmapPh_file, fullfile(workdir, fileBaseName(FmapPh_file)));

FmapMag.img = FmapMag.img(:,:,:,1); % only keep the first volume
FmapMag.hdr.dime.dim(1) = 3;
FmapMag.hdr.dime.dim(5) = 1;

% Convert fieldmap from radians/sec to Hz (fmap/2*pi)
FniiPhHz = FmapPh;
FniiPhHz.img = double(FniiPhHz.img)./(2*pi);
clear FmapPh

switch lower(PED)
   case 'x'
      distort_dim = 1;
      
   case 'x-'
      distort_dim = 1;
      FniiPhHz.img = -1 * FniiPhHz.img;
      
   case 'y'
      distort_dim = 2;
      
   case 'y-'
      distort_dim = 2;
      FniiPhHz.img = -1 * FniiPhHz.img;
      
   otherwise
      error('Unknown PED input: %s', PED)
end


% Distort the fmap magnitude image
fprintf('\nDistorting Fieldmap Magnitude slice: ');
for nslice = 1:size(FniiPhHz.img,3),
   fprintf(', %d', nslice);
   FmapMag.img(:,:,nslice) = EPI_distort_fieldmap_image(FmapMag.img(:,:,nslice), FniiPhHz.img(:,:,nslice), echo_space, distort_dim, 32);
end

% Save distorted magnitude file,
fprintf('\nSaving outputs...');
save_untouch_nii_gz(FmapMag, [output_filebase '.fieldmap_mag.distorted']);


% load dwi, throw away weighted volumes, normalize intensity
dwi = reorient_nifti_sform(dwi_file, fullfile(workdir, fileBaseName(dwi_file)));
enc_out = checkDiffusionEncodingScheme(bvec, bval);
dwi_mask = maskDWI(dwi, fullfile(workdir, 'dwi_mask'), enc_out.bvec, enc_out.bval);
dwi_mask = load_untouch_nii_gz(dwi_mask);
dwi.img = dwi.img(:,:,:,enc_out.zero_bval_mask);
dwi.img = normalize_intensity(dwi.img(:,:,:,1), [1 98], dwi_mask.img);
dwi.hdr.dime.dim(1) = 3;
dwi.hdr.dime.dim(5) = 1;


% interp to dwi coordinates
FmapMag_Dcoord = interp3_nii(FmapMag, dwi, [output_filebase '.fieldmap_mag.distorted.D_coord']);

% normalize fmap to [0 1]
FmapMag_mask = mask_head_pseudo(FmapMag_Dcoord);
FmapMag_Dcoord.img = normalize_intensity(FmapMag_Dcoord.img, [1 98], FmapMag_mask.img);

% write png
overlay_volumes2png(FmapMag_Dcoord.img, dwi.img, [0 1], [output_filebase '.png']);
overlay_volumes2png(FmapMag_Dcoord.img, dwi.img, [0 1], [output_filebase '_edge.png'], 'edge');
rmdir(workdir,'s');
fprintf('\nOverlay successfully saved\n');

end

