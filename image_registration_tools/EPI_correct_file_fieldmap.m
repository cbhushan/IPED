% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function EPI_correct_file_fieldmap(epi_filename, fieldmap_filename, epi_output_filename, options)
% Fieldmap-based distortion correction of EPI images using pixel-shift or Least Square formulation.
% This function is intended to work with EPI data acquired with single (and same) phase encoding
% direction for all EPI slices. Fieldmap should be in units of radians/sec (or one saved by siemens
% scanner).   
%
% Most options are mandatory, see below.

defaultoptions = struct(...
   'phase_encode_direction', 'y', ... x/y/z/x-/y-/z- Direction of distortion
   'echo_spacing', [], ... in Sec
   'siemens_correct', false, ... set to true when the fieldmap file is saved by Siemens (Trio) scanner
   'smooth', false, ... If to apply 3D smoothing to fieldmap before correction
   'smooth_sigma', 0.7, ... in mm. typicaly in range of 0-3.0
   'intensity_correct', true, ... NOTE: If true, noise may not be uniformly distributed in output
   'interpolation_method', 'linear', ... linear / cubic / nearest
   'b0_file', [], ...
   'mask', [], ...
   'checkFOV', true, ... when false, completely skips FOV checks
   'ignore_FOV_errors', false, ... when true, proceeds with correction even when FOV do not match
   'leastsq_sol', true, ... Use Least square formulation 
   'lsqr_scale', 32, ... scaling factor for 
   'lsqr_alpha', 0.27 ... weight for spatial regularization
   );

tags = fieldnames(defaultoptions);
for i=1:length(tags)
   if(~isfield(options,tags{i})),  options.(tags{i})=defaultoptions.(tags{i}); end
end

if length(fieldnames(defaultoptions)) < length(fieldnames(options))
   warning('BDP:UnknownOptions', 'Unknown or too many options found in input option structure');
end

workDir = tempname;
mkdir(workDir);

epi_in = load_untouch_nii_gz(epi_filename, true, workDir);
fmap_in = load_untouch_nii_gz(fieldmap_filename, true, workDir);

% check overlap of fieldmap and epi_data
fm_scanner = checkScanFOV(epi_in, fmap_in, epi_output_filename, options);

switch options.phase_encode_direction
   case {'x', 'x-'}
      permute_vec = [1 2 3 4];
   case {'y', 'y-'}
      permute_vec = [2 1 3 4];
   case {'z', 'z-'}
      permute_vec = [3 2 1 4];
   otherwise
      error('Unknown options.phase_encode_direction')
end

% Siemens correction; fmap_scanner is data after applying nifti header scaling (scl_slope & scl_inter)
% fmap_hz = fmap_scanner/20; 
if options.siemens_correct   
   fm_rad = (double(fm_scanner)./20)*(2*pi);
else
   fm_rad = double(fm_scanner);
end
clear fm_scanner

if options.smooth
   options.smooth_sigma = options.smooth_sigma/min(epi_in.hdr.dime.pixdim(2:4));
   fm_rad = imgaussian(double(fm_rad), options.smooth_sigma);
end

% save the fieldmap which will be used for correction 
fm_rad_nii = epi_in;
fm_rad_nii.hdr.dime.scl_slope = 0;
fm_rad_nii.hdr.dime.scl_inter = 0;
fm_rad_nii.img = fm_rad;
fm_rad_nii.hdr.dime.dim(1) = 3;
fm_rad_nii.hdr.dime.dim(5) = 1;
save_untouch_nii_gz(fm_rad_nii, fullfile(fileparts(epi_output_filename), 'fieldmap.radians.used.nii.gz'), 16, workDir);
clear fm_rad_nii

% permute data such that 1st dimension is phase encode direction
epi_in.img = permute(epi_in.img, permute_vec);
fm_rad = permute(fm_rad, permute_vec);
epi_res = permute(epi_in.hdr.dime.pixdim(2:5), permute_vec);

switch options.phase_encode_direction
   case {'x-', 'y-', 'z-'}
      fm_rad = -1*fm_rad;
end

% correct EPI images using pixelshift
if options.leastsq_sol
   msg = {'\n', 'Computing initial estimate for distortion correction using Least Square formulation...\n'};
else
   msg = {'\n', 'Correcting EPI (diffusion) image for distortion using pixel-shift...\n'};
end
fprintf(bdp_linewrap(msg));

sz = size(epi_in.img);
if length(sz)<4
   sz(4) = 1;
end

epi_out = epi_in;
epi_out.img = epi_out.img.*0;

cpb = ConsoleProgressBar(); % Set progress bar parameters
cpb.setMinimum(1); cpb.setMaximum(sz(4)); cpb.start();
for k = 1:sz(4)
   epi_out.img(:,:,:,k) = EPI_correct_fieldmap_pixelshift(epi_in.img(:,:,:,k), epi_res, fm_rad/(2*pi), ...
      options.echo_spacing, options.intensity_correct, options.interpolation_method);
   
   msg = sprintf('%d/%d volumes done', k, sz(4));
   cpb.setValue(k); cpb.setText(msg);
end
cpb.stop();


% correct EPI images using LSQR
if options.leastsq_sol
   msg = {'\n', ['Correcting EPI (diffusion) images for distortion using Least Square formulation...']};
   fprintf(bdp_linewrap(msg));
   
   epi_lsqr_corr = permute(epi_in.img,[1 2 4 3]).*0;
   cpb = ConsoleProgressBar(); % Set progress bar parameters
   cpb.setMinimum(1); cpb.setMaximum(sz(3)); cpb.start();
   phase_dir_ind = ones(sz(4),1);
   
   for parslc = 1:sz(3)
      data_slc = double(squeeze(epi_in.img(:,:,parslc,:)));
      deltaB0_slice = squeeze(fm_rad(:,:,parslc))/(2*pi);
      x_init = abs(squeeze(epi_out.img(:,:,parslc,:)));
      epi_lsqr_corr(:,:,:,parslc) = EPI_correct_slice_invLeastSq_image(data_slc, phase_dir_ind, deltaB0_slice, ...
         options.echo_spacing, options.lsqr_scale, options.lsqr_alpha, false, sz(4), x_init);
      msg = sprintf('%d/%d slices done', parslc, sz(3));
      cpb.setValue(parslc); cpb.setText(msg);
   end
   cpb.stop();
   epi_out.img = permute(epi_lsqr_corr, [1 2 4 3]);
   clear epi_lsqr_corr data_slc deltaB0_slice x_init
end


% back to original size
epi_out.img(epi_out.img<0) = 0;
epi_out.img = permute(epi_out.img, permute_vec);

fprintf('\nSaving file...'); drawnow('update');
epi_out.hdr.dime.scl_slope = 0;
epi_out.hdr.dime.scl_inter = 0;
save_untouch_nii_gz(epi_out, epi_output_filename, workDir);
fprintf('Done')

rmdir(workDir, 's');
end


function fmap_out = checkScanFOV(epi_in, fmap_in, epi_output_filename, options)
% check overlap of fieldmap and epi_data and returns the fieldmap data resampled in coordinate of
% epi_data

vox_diff_allowed = 1; % voxel difference allowed in FOV
[~, X_vol1, Y_vol1, Z_vol1, res1, Tvol1] = get_original_grid_data(epi_in);
[~, ~, ~, ~, res2, Tvol2] = get_original_grid_data(fmap_in);

if options.checkFOV && (~isequal(size(epi_in.img), size(fmap_in.img)) || norm(Tvol1(:)-Tvol2(:))>1e-2)
   
   overlay_fname = [remove_extension(epi_output_filename) '.fielmap_overlay'];
   msg = {'\n', ['Fieldmap and EPI (diffusion) scan have different headers/acquisition scheme. '...
      'PNG images showing their overlay will be generated with name: ' escape_filename(overlay_fname)], ...
      '\n', ['It is *highly* recommended to check the overlay images to make sure fieldmap and EPI data '....
      'overlap correctly. Input fieldmap scan should be pre-registered to EPI (diffusion) image for '...
      'fieldmap-based correction. '], '\n'};
   fprintf(bdp_linewrap(msg));
   
   if ~isempty(options.mask)
      msk = load_untouch_nii_gz(options.mask);
   else
      sz = size(epi_in.img);
      msk.img = true(sz(1:3));
   end
   
   if ~isempty(options.b0_file)
      b0 = load_untouch_nii_gz(options.b0_file, true);
   else
      b0.img = epi_in.img(:,:,:,1);
   end
   
   c(1,:) = X_vol1(msk.img>0);
   c(2,:) = Y_vol1(msk.img>0);
   c(3,:) = Z_vol1(msk.img>0);
   c(4,:) = 1;
   c = Tvol2\c;
   c_min = transpose(min(c,[],2)); c_min(end) = [];
   c_max = transpose(max(c,[],2)); c_max(end) = [];
   
   fmap_out = myreslice_nii(fmap_in, 'linear', X_vol1, Y_vol1, Z_vol1);
   fmap_sz = size(fmap_in.img);
   
   fm_norm = normalize_intensity(fmap_out, [8 92], msk.img>0);
   epi_norm = normalize_intensity(b0.img(:,:,:,1), [8 96], msk.img>0);
   overlay_volumes2png(fm_norm, epi_norm, [0 1], overlay_fname, 'rview', [85 100]);
   clear fm_norm epi_norm
   
   if min(c_min)<(-vox_diff_allowed) || min(fmap_sz-c_max)<(1-vox_diff_allowed) % Allow one voxel differences due to partial volume errors
      num_vox_outside = sum( c(1,:)<-0.5 | c(2,:)<-0.5 | c(3,:)<-0.5 |...
         c(1,:)>fmap_sz(1)-0.5 | c(2,:)>fmap_sz(2)-0.5 | c(3,:)>fmap_sz(3)-0.5 );
      
      if options.ignore_FOV_errors
         msg = {'\n', ['EPI (diffusion) field of view (FOV) is not totally covered by input fieldmap. '...
            'Number of voxels in EPI (diffusion) scans outside FOV of fieldmap: ' num2str(num_vox_outside) '. ' ...
            'But --ignore-fieldmap-FOV is detected. So, BDP will ignore it and continue. '...
            'However, BDP has saved some images (.png files) which show the overlay of EPI (diffusion) data and fieldmap. '...
            'It is *highly* recommended to check the overlay images and make sure fieldmap and EPI (diffusion) data '...
            'overlap resonably. Input fieldmap scan should be pre-registered to EPI image for fieldmap-'...
            'based correction.\n']};
         fprintf(bdp_linewrap(msg));
      else
         msg = ['\nEPI (diffusion) field of view (FOV) is not totally covered by input fieldmap. '...
            'Number of voxels in EPI (diffusion) scans outside FOV of fieldmap: ' num2str(num_vox_outside) '. ' ...
            'BDP has saved some images (.png files) which show the overlay of EPI (diffusion) data and fieldmap. '...
            'Please check these overlay images and make sure that fieldmap file has correct header. '...
            'Input fieldmap scan should be pre-registered to EPI (diffusion) image for fieldmap-based correction. '...
            'If you think that you got this error in mistake then you can suppress this error by '...
            're-running BDP and appending flag --ignore-fieldmap-FOV.\n'];
         error('BDP:InsufficientFieldmapFOV', bdp_linewrap(msg));
      end
   end
   
else
   fmap_out = myreslice_nii(fmap_in, 'linear', X_vol1, Y_vol1, Z_vol1);
end
end


function [epi_out, shift] = EPI_correct_fieldmap_pixelshift(epi_vol, epi_res, deltaB0, echo_space, intensity_correct, method)
% epi_vol - 3D EPI volume, with 1st dimension as phase encode & second dimension as readout direction
% deltaB0 - in Hz; 3D volume (dimensions same as epi_vol)
% echo_spacing  - in sec
% epi_res - vector of length three

sz = size(epi_vol);
shift = double(deltaB0) * sz(1) * double(echo_space); % in voxels

% warp EPI
[x_g, y_g, z_g] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
x_warp = x_g + shift;
epi_out = interpn(double(epi_vol), x_warp, y_g, z_g, method, 0);

if intensity_correct
   [~, grad_x, ~] = gradient(shift*epi_res(1), epi_res(1), epi_res(2), epi_res(3));
   epi_out = epi_out.*(1+grad_x);
end
end
