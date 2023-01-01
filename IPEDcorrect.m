% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [outfile, dwi_corr] = IPEDcorrect(IPED_config_filename, fieldmap_filename, echo_spacing, output_base, user_opts)
% Estimates distortion corrected diffusion weighted images from multiple phase encoded data.
% Either of following syntax is supported:
% 
%   IPEDcorrect(IPED_config_filename, fieldmap_filename, echo_spacing, output_base)
%   IPEDcorrect(IPED_config_filename, fieldmap_filename, echo_spacing, output_base, opts)
% 
% `IPED_config_filename` is a string specifying the file-name of the IPED-config file. See below for more details about format of the IPED-config file.
% 
% `fieldmap_filename` is a string specifying the file-name of the fieldmap saved in NIfTI-1 format. The fieldmap can acquired or estimated using any method but it must be specified in units of radians/sec. The voxel resolution of fieldmap need not match to that of the diffusion images, however for best correction it should cover the field-of-view of the diffusion scans.
% `echo_spacing` is a number specifying the echo-spacing in units of secs. Eg: For an echo spacing of 0.36ms, use 0.00036.
% `output_base` is a string specifying the file-prefix for the outputs written by IPEDcorrect. Eg: when output_base is set to sub02.IPED-corrected, the file-names of the outputs written by IPEDcorrect would be sub02.IPED-corrected.nii.gz, sub02.IPED-corrected.bvec and likewise.
% `opts` is an optional input which specifies several optional arguments, like regularization parameters, for distortion correction. 
%
% Please see README and USAGE files for more details.
%
% Current Limitations: 
%   * Each PED can appear only once in conf-file.
%   * No automatic detection or check/test for full-data
%   * Spatial regularization is not aware of spatial resolution (OK for isotropic resolution)
%

print_IPED_preamble();

if nargin<4 || nargin>5
   link = 'https://github.com/cbhushan/IPED';
   msg = sprintf('\n%s\n\t %s', bdp_linewrap(['Incorrect number of input arguments. IPEDcorrect() needs either '...
      '4 or 5 input arguments. Please see usage details at:']), link);
   error(msg);
end

opts = struct(...
   'spatial_beta', 0.27, ...
   'spherical_alpha', 0.04, ...
   'shOrder', 8, ... must be even >0
   'sph_wt', 'G_mag', ... 'G_mag'/[]/filename
   'G_mag_low', 1*2*pi, ...rad/sec/mm
   'G_mag_high', 9*2*pi, ...rad/sec/mm
   'estimation_mask', '', ... mask filename
   'isfull', false, ...
   'deform_upscale', 16, ...
   'save_SH_coeff', 1 ... 0/1/2; 1=BrainSuite style SH-coeff; 2=single 4D file
   );

if exist('user_opts','var')
   tags = fieldnames(opts);
   for i=1:length(tags)
      if isfield(user_opts, tags{i})
         opts.(tags{i}) = user_opts.(tags{i});
      end
   end
   if any(~ismember(fieldnames(user_opts), tags))
      warning('BDP:UnknownOptions', 'unknown options found in user_opts');
   end
end

workdir = tempname();
mkdir(workdir);
output_base =  remove_extension(output_base);


% read IPED-text file
fprintf('\nReading input files...')
s = IPEDtxtReader(IPED_config_filename);

if length(s)==1   
   msg = ['Insufficient phase encodings! Atleast two different phase encodings are required '...
      'for interlaced correction. Only one dataset was read from IPED-conf file.'];
   error(bdp_linewrap(msg));
end

% Find orientation of slices - Also see PED_dir() below
PEDs_avl = {1, length(s)};
PEDs_avl_non_polar = {1, length(s)};
for k = 1:length(s)
   PEDs_avl_non_polar{k} = strrep(s(k).PED, '-' ,'');
   PEDs_avl{k} = s(k).PED;
end
PEDs_avl_non_polar = unique(PEDs_avl_non_polar); % unique PEDs ignoring polarity
PEDs_avl = unique(PEDs_avl);

if length(PEDs_avl)==1
   msg = ['Insufficient phase encodings! Atleast two different phase encodings are required '...
      'for interlaced correction. Only one PED direction was detected in IPED-conf file.'];
   error(bdp_linewrap(msg));
   
elseif length(PEDs_avl)~=length(s)
   msg = ['Unsupported format of input files! Each unique phase encoding direction can be '...
      'specified only once in IPED-conf file - Please concatenate data for each phase-encoding' ...
      'into one set of nifti/bvec/bval/bmat file and try again.'];
   error(bdp_linewrap(msg));

elseif length(PEDs_avl_non_polar)>2
   msg = ['Unsupported acquisition - EPI plane for all DWIs acquired with different PED '...
      'must coincide (i.e. ignoring ' ...
      'polarity of PED, there could be maximum of two PEDs).'];
   error(bdp_linewrap(msg));
end

possible_PEDs = {'x' 'y' 'z'};
PEDs_ind = find(ismember(possible_PEDs, PEDs_avl_non_polar));
PEDs_avl_non_polar = possible_PEDs(PEDs_ind); % unique PEDs ignoring polarity - in sorted order
temp = setdiff([1 2 3], PEDs_ind, 'stable');
permute_vec = [PEDs_ind temp 4]; % permutation order of vols (after rotating to RAS)



% viz directions
cout = randCMapHue(length(s));
h = displayDiffusionDirections(s(1).bvec, cout(1,:));
for k = 2:length(s)
   h = displayDiffusionDirections(s(k).bvec, cout(k,:), h);
end
fname = [output_base '.bvecs.png']; 
drawnow(); state = pause('query'); pause on; pause(0.5); pause(state); % some dirty hack!
saveas(h, fname, 'png');
fprintf('\nSaved image showing diffusion encoding directions: %s\n', fname);


% reorient to RAS, sanity checks on files
fprintf('\nChecking input nifti files...')
[~, tempfile] = check_nifti_file(s(1).dwi_file, workdir);
outFile = fullfile(workdir, [fileBaseName(s(1).dwi_file) '.RSA.nii.gz']);
[~, ~, ~, nii] = reorientDWI_BDP(tempfile, outFile, s(1).bmat, s(1));
sz = size(nii.img);
res = abs(nii.hdr.dime.pixdim(2:4));

if length(sz)~=4
   error(bdp_linewrap(['All the input images must be 4D-nifti volumes. 4th dimension '...
      'should represent the diffusion-weighting dimension.\n']))
end

for k = 1:length(s)
   [~, tempfile] = check_nifti_file(s(k).dwi_file, workdir);
   outFile = fullfile(workdir, [fileBaseName(s(k).dwi_file) '.RSA.nii.gz']);
   [~, ~, ~, nii, bmat, bvec] = reorientDWI_BDP(tempfile, outFile, s(k).bmat, s(k));
   
   sz_tmp = size(nii.img);
   if ~isequal(sz_tmp(1:3), sz(1:3))
      error('\nSize of the image is not consistent with others: %s\n', s(k).dwi_file);
   end
   
   if ~isempty(bvec) && ~isequal(sz_tmp(4), size(bvec, 1))
      error('\nNumber of DWIs is not consistent with specified diffusion (bvec) parameter: %s\n', s(k).dwi_file);
   end
   
   if ~isempty(bmat) && ~isequal(sz_tmp(4), size(bmat, 3))
      error('\nNumber of DWIs is not consistent with specified diffusion (bmat) parameter: %s\n', s(k).dwi_file);
   end
   
   % add RAS data
   s(k).dwi_RAS = nii;
   s(k).bvec_RAS = bvec;
   s(k).bmat_RAS = bmat;
end

% reorient fieldmap to RAS
[~, tempfile] = check_nifti_file(fieldmap_filename, workdir);
fieldmap_filename = fullfile(workdir, [fileBaseName(fieldmap_filename) '.RSA.nii.gz']);
reorient_nifti_sform(tempfile, fieldmap_filename);



% correct b=0 using single PED - for user's check of PED direction
fprintf('\n\nSingle-PED correction - For checking PED in config files...\n')
b0_corr_1PED = '';
for k = 1:length(s)
   nii = s(k).dwi_RAS;
   nii.img = nii.img(:,:,:,s(k).zero_bval_mask);
   nii.hdr.dime.dim(5) = size(nii.img, 4);
   tempfile = save_untouch_nii_wrapper(nii, fullfile(workdir, [fileBaseName(s(k).dwi_file) '.RSA.0_diffusion.nii.gz']));
   epi_output_filename = [output_base '.' s(k).PED '.singlePED.0_diffusion.correct.nii.gz'];
   
   fmap_opts.phase_encode_direction = s(k).PED;
   fmap_opts.echo_spacing = echo_spacing;
   fmap_opts.leastsq_sol = false;
   fmap_opts.ignore_FOV_errors = true;
   EPI_correct_file_fieldmap(tempfile, fieldmap_filename, epi_output_filename, fmap_opts);
   delete(tempfile);
   b0_corr_1PED = [b0_corr_1PED '\t' epi_output_filename '\n'];
end
msg = {'\n\n', ['Please check if the following distortion corrected (b=0) images look reasonable. '...
   'These are corrected using a single image with one-PED (& are not final output). If these images '...
   'do not look appropriately corrected for distortion then the final output would most likely be '...
   'inaccurately corrected. In such situation, the input PED direction in IPED-conf '...
   'file may need to be setup correctly - in most cases changing the polarity of PED should do the trick '...
   '(eg. replace x by x- and vice-versa).\n\n']};
fprintf(bdp_linewrap(msg));
fprintf(b0_corr_1PED);
clear nii


% Get fieldmap in EPI-coord
target_file = [output_base '.' s(1).PED '.singlePED.0_diffusion.correct.nii.gz'];
file_out = [output_base '.fieldmap.radians.used.nii.gz'];
fieldmap_filename = fullfile(fileparts(output_base), 'fieldmap.radians.used.nii.gz'); % Used fieldmap, saved by EPI_correct_file_fieldmap()
fmap_rad = interp3_nii(fieldmap_filename, target_file, file_out, 'linear');
fmap_hz = double(fmap_rad.img)/(2*pi);
fmap_hz = permute(fmap_hz, permute_vec);
delete(fieldmap_filename); % It was saved by EPI_correct_file_fieldmap()
clear fmap_rad


% reorient masks and weights to RAS
if isempty(opts.estimation_mask)
   est_mask = true(sz(1:3));
else
   hdr_src_file = [output_base '.' s(1).PED '.singlePED.0_diffusion.correct.nii.gz'];
   output_file = fullfile(workdir, 'estimation_mask.nii.gz');
   est_mask = fixBSheader(hdr_src_file, opts.estimation_mask, output_file);
   est_mask = est_mask.img>0;
end

if ~isempty(opts.sph_wt)
   if ischar(opts.sph_wt)
      if ~strcmpi(opts.sph_wt, 'G_mag')
         hdr_src_file = [output_base '.' s(1).PED '.singlePED.0_diffusion.correct.nii.gz'];
         output_file = fullfile(workdir, 'spherical_weight.nii.gz');
         wt = fixBSheader(hdr_src_file, opts.sph_wt, output_file);
         opts.sph_wt = double(wt.img);
         clear wt         
      end
   else
      error('opts.sph_wt must be either empty or "G_mag" or filename of the weight-file.\n')
   end
end


% permute dims for distortion correction
for k = 1:length(s)
   s(k).dwi_permute = permute(s(k).dwi_RAS.img, permute_vec);
   s(k).dwi_RAS.img = []; % clear to save memory
   
   s(k).b0_permute = mean(s(k).dwi_permute(:,:,:,s(k).zero_bval_mask), 4);
   s(k).dwi_permute = s(k).dwi_permute(:,:,:, ~s(k).zero_bval_mask);
   
   temp = load_untouch_nii_gz([output_base '.' s(k).PED '.singlePED.0_diffusion.correct.nii.gz']);
   temp = permute(temp.img, permute_vec);
   s(k).b0_init_permute = mean(temp, 4);   
end
est_mask = permute(est_mask, permute_vec);
if isnumeric(opts.sph_wt),
    opts.sph_wt = permute(opts.sph_wt, permute_vec);
end

% correct b=0 using all data
fprintf('\n\nCorrecting b=0 using all PED image...\n')
sz_perm = sz(permute_vec);
res_perm = res(permute_vec(1:3)); 
b0_full_corr = zeros(sz_perm(1:3));
isfull = true;
num_diffusion_dir = 1;

cpb = ConsoleProgressBar(); % Set progress bar parameters
cpb.setMinimum(1); cpb.setMaximum(sz_perm(3)); cpb.start();
for slc = 1:sz_perm(3)   
   data_slc = zeros([sz_perm(1:2) length(s)]);
   b0_init = zeros(sz_perm(1:2));
   phase_dir = zeros(length(s), 1);
   for k = 1:length(s)
      data_slc(:,:,k) = double(s(k).b0_permute(:,:,slc));
      b0_init = b0_init + double(s(k).b0_init_permute(:,:,slc));
      phase_dir(k) = PED_dir(s(k).PED, PEDs_avl_non_polar);
   end
   b0_init = b0_init/length(s);
   deltaB0_slice = squeeze(fmap_hz(:,:,slc));
   
   b0_full_corr(:,:,slc) = EPI_correct_slice_invLeastSq_image(data_slc, phase_dir, deltaB0_slice, ...
      echo_spacing, opts.deform_upscale, opts.spatial_beta, isfull, num_diffusion_dir, ...
      b0_init); %, mask_undistorted, mask_dwi_slice)
   
   msg = sprintf('%d/%d slices done', slc, sz_perm(3));
   cpb.setValue(slc); cpb.setText(msg);
end
cpb.stop();
dwi_corr.b0_corr = ipermute(b0_full_corr, permute_vec);


% initialization for IPED correction using pixel-shift



% IPED correction
fprintf('\n\nCorrecting DWIs using IPED...\n')
sz_perm = sz(permute_vec);
opts.G_mag_low = opts.G_mag_low/2/pi; % in Hx/mm
opts.G_mag_high = opts.G_mag_high/2/pi;
reg_param = [opts.spherical_alpha opts.spatial_beta];

bvecs = [];
for k = 1:length(s)
   bvecs = cat(1, bvecs, s(k).bvec_RAS(~s(k).zero_bval_mask, :));
end

if opts.isfull
   num_diffusion_dir = sum(~s(1).zero_bval_mask);
else
   num_diffusion_dir = size(bvecs, 1);
end
d_corr = zeros([sz_perm(1:3) size(bvecs, 1)]);
sh_corr = zeros([sz_perm(1:3) (opts.shOrder+1)*(opts.shOrder+2)/2]);

cpb = ConsoleProgressBar(); % Set progress bar parameters
cpb.setMinimum(1); cpb.setMaximum(sz_perm(3)); cpb.start();
for slc = 1:sz_perm(3)
   dwi_slc = zeros([sz_perm(1:2) size(bvecs, 1)]);
   phase_dir = zeros(size(bvecs, 1), 1);
   c = 1;
   for k = 1:length(s)
      n = sum(~s(k).zero_bval_mask);
      dwi_slc(:,:, c:(c+n-1)) = double(squeeze(s(k).dwi_permute(:,:,slc,:)));
      phase_dir(c:(c+n-1)) = PED_dir(s(k).PED, PEDs_avl_non_polar);
      c = c+n;
   end
   deltaB0_slice = squeeze(fmap_hz(:,:,slc));
   
   opts_IPED.spherical_weight = [];
   opts_IPED.mask_undistorted = est_mask(:,:,slc);
   if ~isempty(opts.sph_wt)
      if ischar(opts.sph_wt) && strcmpi(opts.sph_wt, 'G_mag')
         [gy, gx] = gradient(deltaB0_slice, res_perm(1), res_perm(2));
         wt = sqrt(gx.^2 + gy.^2);
         wt(wt<opts.G_mag_low) = opts.G_mag_low;
         wt(wt>opts.G_mag_high) = opts.G_mag_high;
         opts_IPED.spherical_weight = wt(:);
         clear wt
      elseif isnumeric(opts.sph_wt)
         opts_IPED.spherical_weight = vect(opts.sph_wt(:,:,slc));
      end
   end
   
   [d_corr(:,:,slc,:), sh_corr(:,:,slc,:)] = EPI_correct_slice_invLeastSq_spharm(dwi_slc, phase_dir, deltaB0_slice, echo_spacing, ...
      opts.deform_upscale, bvecs, opts.shOrder, reg_param, opts.isfull, num_diffusion_dir, opts_IPED);
   
   msg = sprintf('%d/%d slices done', slc, sz_perm(3));
   cpb.setValue(slc); cpb.setText(msg);
end
dwi_corr.d_corr = ipermute(d_corr, permute_vec);
dwi_corr.sh_corr = ipermute(sh_corr, permute_vec);
dwi_corr.nii_struct = s(1).dwi_RAS;


% save corrected outputs
fprintf('\n\nSaveing output files...')
nii = s(1).dwi_RAS;
nii.img = cat(4, dwi_corr.b0_corr, dwi_corr.d_corr);
nii.hdr.dime.dim(1) = 4;
nii.hdr.dime.dim(5) = size(nii.img, 4);
outfile.dwi_file = save_untouch_nii_gz(nii, [output_base '.RAS.IPED_correct.nii.gz'], 16, workdir);

bvec_file = [output_base '.RAS.IPED_correct.bvec'];
bval_file = [output_base '.RAS.IPED_correct.bval'];
bmat_file = [output_base '.RAS.IPED_correct.bmat'];
if opts.isfull
   writeBvecBvalFile(s(1).bvec_RAS, bvec_file, s(1).bval, bval_file);
   writeBmatFile(s(1).bmat_RAS, bmat_file);
else
   bvals = [];
   bmats = [];
   for k = 1:length(s)
      bvals = cat(1, bvals, s(k).bval(~s(k).zero_bval_mask, :));
      bmats = cat(3, bmats, s(k).bmat_RAS(:,:, ~s(k).zero_bval_mask));
   end
   bvals = [0; bvals];
   bvecs = cat(1, [0 0 0], bvecs);
   bmats = cat(3, zeros(3), bmats);
   
   writeBvecBvalFile(bvecs, bvec_file, bvals, bval_file);
   writeBmatFile(bmats, bmat_file);
end
outfile.bvec_file = bvec_file;
outfile.bval_file = bval_file;
outfile.bmat_file = bmat_file;


% write SH coefficients
if opts.save_SH_coeff==1
   pdir = fullfile(fileparts(output_base), 'DWI_IPED_SHcoeff');
   mkdir(pdir);
   outfile.SHc_dir =  pdir;
   
   SH_file_base = fullfile(pdir, [fileBaseName(output_base) '.RAS.IPED_correct.SH']);
   odf_file = [SH_file_base '.odf']; % BrainSuite SHc description file
   outfile.SHc_description_file = odf_file;
   odf_fid = fopen(odf_file, 'w');
   
   nii = s(1).dwi_RAS;
   nii.hdr.dime.dim(1) = 3;
   nii.hdr.dime.dim(5) = 1;
   for k = 1:size(dwi_corr.sh_corr, 4)
      fname = sprintf('%s.%03d.nii.gz', SH_file_base, k);
      [~, name, ext] = fileparts(fname);
      fprintf(odf_fid, '%s\n', [name, ext]);
      nii.img = dwi_corr.sh_corr(:,:,:,k);
      save_untouch_nii_gz(nii, fname, 16);
   end
   
elseif opts.save_SH_coeff==2 % 4D SH-coeff file
   nii = s(1).dwi_RAS;
   nii.img = dwi_corr.sh_corr;
   nii.hdr.dime.dim(1) = 4;
   nii.hdr.dime.dim(5) = size(nii.img, 4);
   outfile.SHc_4D_file = save_untouch_nii_gz(nii, [output_base '.IPED_correct.nii.gz'], 16);
end
fprintf('Done.\n')

rmdir(workdir, 's');
end



function d = PED_dir(PED, PEDs_avl_non_polar)
% Returns an integer which specified 'phase_dir' for invLeastSq functions
%   PEDs_avl - cell array of unique PEDs ignoring polarity
%   PED - PED string

if length(PEDs_avl_non_polar)==1
   % Only one PED - the fist dim of the volume after permuting would be the PED
   if ismember(PED, {'x', 'y', 'z'})
      d = 1;
   else
      d = 2;
   end
   
elseif length(PEDs_avl_non_polar)==2
   % Two PEDs - check position in PEDs_avl
   ind = find(ismember(PEDs_avl_non_polar, strrep(PED, '-' ,'')));
   
   if ind==1 % fist dim after permuting
      if ismember(PED, {'x', 'y', 'z'})
         d = 1;
      else
         d = 2;
      end
      
   else     % second dim after permuting
      if ismember(PED, {'x', 'y', 'z'})
         d = 3;
      else
         d = 4;
      end
   end
else
   error('This is not supported. You should not be here!')
end

end
