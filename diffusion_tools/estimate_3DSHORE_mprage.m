% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function estimate_3DSHORE_mprage(data_file, target_file, affinematrixfile, bMatrices, user_options)
% Computes spherical harmonic coefficients of ODFs using 3DSHORE. Saves the output in
% BrainSuite friendly format. 
%
%  No eddy current or motion correction.
%
%  data_file - 4D nifti diffusion file name
%  target_file - bfc file name (this will be also used as mask)
%  affinematrixfile - file name of .mat file saved after rigid registration of diffusion and
%                     bfc file
%  bMatrices - .bmat file name or array of matrices (3x3xN)
%  user_options - optional input for Spherical harmonic parameters (See below for default options)

workdir = tempname;
mkdir(workdir);

% setting up the defaults options
opt = struct( ...
   'RadialOrder', 8, ...
   'shore_lambdaN', 1e-8, ...
   'shore_lambdaL', 1e-8, ...
   'shore_zeta', 700, ...
    'estimate_odf_3DSHORE', true, ...
   'diffusion_time', 2.5*10^-3, ...
   'diffusion_coeff',0.7*10^-3,...
   'SHORE_out_dir', fullfile(fileparts(data_file), '3DSHORE'), ...
   'mprage_coord_suffix', '.T1_coord', ...
   'bval_ratio_threshold', 45, ...
   't1_mask_file', target_file ...
    );

fprintf('\nEstimating 3DSHORE ODFs in T1-coordinate...');

% set options using user_options
if exist('user_options', 'var')
   if isfield(user_options, 'file_base_name')
      opt.SHORE_out_dir = fullfile(fileparts(user_options.file_base_name), '3DSHORE');
   end
   
   fnames = fieldnames(opt);
   for iname = 1:length(fnames)
      if isfield(user_options, fnames{iname})
         opt = setfield(opt, fnames{iname}, getfield(user_options, fnames{iname}));
      end
   end
end


fname = fileBaseName(data_file);
SHORE_output_file_base = fullfile(opt.SHORE_out_dir, fname);
if opt.estimate_odf_3DSHORE && ~exist(opt.SHORE_out_dir, 'dir'), mkdir(opt.SHORE_out_dir); end

% Apply the rigid transform to header of dwi file
load(affinematrixfile);
data_file_transformed = fullfile(workdir, [Random_String(15) '.nii.gz']);
[~,Tnew] = affine_transform_nii(data_file, M_world, origin, data_file_transformed);
clear Tnew;

% load bmatrices
if ischar(bMatrices)
   bMatrices = readBmat(bMatrices);
end

% load data and corresponding grid points
[dataIn, ~, ~, ~, res_dwi, Tdwi] = get_original_grid_data(data_file_transformed);
[~, X_target, Y_target, Z_target, res_target, T_target] = get_original_grid_data(target_file);

dwi_vol_size = size(dataIn);
target_vol_size = size(X_target);
nDir = dwi_vol_size(4);

if nDir~=size(bMatrices, 3)
   error('BDP:InconsistentDiffusionParameters', ['Number of diffusion image should be exactly same' ...
      'as number of bMatrices/bvec (in .bmat/.bvec file, if any).']);
end

% Rotate bmatrices
for iDir = 1:nDir 
   Rmat = (T_target*diag(1./[res_target 1])) \ (Tdwi*diag(1./[res_dwi 1]));
   Rmat = Rmat(1:3, 1:3);
   bMatrices(:,:,iDir) = Rmat*bMatrices(:,:,iDir)*(Rmat');
end

% Find diffusion encoding direction (same as largest eigen vector of bMatrices)
DEout = checkDiffusionEncodingScheme(bMatrices, opt.bval_ratio_threshold);
ind = 1:nDir; %find(~DEout.zero_bval_mask);
q = zeros(numel(ind),3);
del_t = opt.diffusion_time;
for i = 1:numel(ind)
    [V,D] = eig(bMatrices(:,:,ind(i)));
    [maxb,tmp] = max(diag(D));
    q(i,:) = V(:,tmp);
    qrad(i) = sqrt(maxb/(4*pi^2*del_t));
end

% template 3D volume
vd = load_untouch_nii_gz(target_file,true);
temp3 = vd;
temp3.img = [];
temp3.hdr.dime.dim(1) = 3;
temp3.hdr.dime.dim(5) = 1;

vd = load_untouch_nii_gz(opt.t1_mask_file);
targetMask = vd.img>0;
clear vd

% dilate mask
% pix = 2;
% [x1,y1,z1] = ndgrid(-pix:pix);
% se2 = (sqrt(x1.^2 + y1.^2 + z1.^2) <=pix);
% targetMask = imdilate(targetMask>0, se2)>0;

% get location in original grid of datafile
c(1,:) = X_target(:); clear X_target;
c(2,:) = Y_target(:); clear Y_target;
c(3,:) = Z_target(:); clear Z_target;
c(4,:) = 1;
c = Tdwi\c;
X_dwi_target = reshape(c(1,:), target_vol_size);
Y_dwi_target = reshape(c(2,:), target_vol_size);
Z_dwi_target = reshape(c(3,:), target_vol_size);
clear c

% Voxel indexing starts from 0
[X_dwi, Y_dwi, Z_dwi] = ndgrid(0:dwi_vol_size(1)-1, 0:dwi_vol_size(2)-1,  0:dwi_vol_size(3)-1);

% 3DSHORE Computation
RadialOrder = opt.RadialOrder;
lambdaN = opt.shore_lambdaN;
lambdaL = opt.shore_lambdaL;
Dcoeff = opt.diffusion_coeff;
zeta = 1/(8*pi^2*Dcoeff*del_t);
% zeta = opt.shore_zeta;

% Regularization
Lshore = l_shore(RadialOrder);
Nshore = n_shore(RadialOrder);

% Shore basis
[Sh,R1,Snew,N1] = shore_basis(qrad,[q(:,1),q(:,2),q(:,3)],RadialOrder,zeta);
% Shore basis pseudo inverse "H+" matrix
Sh_LS = (Sh'*Sh + lambdaN.*Nshore + lambdaL.*Lshore )\Sh';

% Radial integral term
[radInt,N] = shore_odf_factor(zeta,RadialOrder);

shoreMatrixReg = Sh_LS'*diag(radInt);

nVoxels = sum(targetMask(:));
sh_3dshore = zeros([max(N) nVoxels], 'single');
sh_coeff = zeros([max(N) target_vol_size(1)*target_vol_size(2)], 'single'); % one slice

cpb = ConsoleProgressBar(); % Set progress bar parameters
cpb.setMinimum(1);     
cpb.setMaximum(target_vol_size(3)); 
cpb.start();
ivox = 0;
target_ind = [];
for slice = 1:target_vol_size(3)
   slice_z_coord = Z_dwi_target(:,:,slice);
   
   if max(slice_z_coord(:))>0 && min(slice_z_coord(:))<=dwi_vol_size(3)     
      msk = targetMask(:,:,slice);
      dwimages = zeros([sum(msk(:)) nDir]);
      
      % compute only inside mask
      X = X_dwi_target(:,:,slice);
      Y = Y_dwi_target(:,:,slice); 
      Z = Z_dwi_target(:,:,slice);
      
      X = X(msk);
      Y = Y(msk);
      Z = Z(msk);
      
      for iDir = 1:nDir
         dwimages(:,iDir) = interpn(X_dwi, Y_dwi, Z_dwi, double(dataIn(:,:,:,iDir)), X, Y, Z, 'linear', 0);
      end
      dwimages(~isfinite(dwimages))=0;

      % compute SHORE coefficients
      b0mean = mean(dwimages(:, DEout.zero_bval_mask),2);
      dwimages = dwimages(:, ind);      
      dwimages = dwimages./b0mean(:,ones(1,size(dwimages,2))); % normalize by b=0 image      
      dwimages(~isfinite(dwimages))=0;
      dwimages = permute(dwimages,[2,1]); % 1st dim is different diffusion weighting
      
      shore_coeff = single(shoreMatrixReg'*dwimages);

      % Convert to SH coefficients
	  if(size(dwimages,2) ~=0)
		  for num =  1:max(N),
			  sh_3dshore(num, ivox+1:ivox+size(dwimages,2)) = sum(shore_coeff(N==num,:),1);
		  end;
      end;
      msk_ind = find(msk);
      target_ind = [target_ind; msk_ind+(slice-1)*prod(target_vol_size(1:2))];
      ivox = ivox + size(dwimages,2);
   end
   
   text = sprintf('%d/%d slices done', slice, target_vol_size(3));
   cpb.setValue(slice); cpb.setText(text); 
end
cpb.stop();
fprintf('\n');
clear dwimages b0mean

shore_fid = fopen([SHORE_output_file_base '.SH.3DSHORE' opt.mprage_coord_suffix '.odf'], 'w');

disp('Writing 3DSHORE ODF files to disk...')
temp3.img = zeros(target_vol_size);
for k = 1:max(N),
    fname = sprintf('%s.SH.3DSHORE.%03d.nii.gz', SHORE_output_file_base, k);
    [~, name, ext] = fileparts(fname);
    fprintf(shore_fid, '%s\n', [name, ext]);
    temp3.img(target_ind) = sh_3dshore(k,1:length(target_ind));
    save_untouch_nii_gz(temp3, fname, 16);
end
fclose(shore_fid);

rmdir(workdir, 's');
fprintf('Estimated 3DSHORE ODFs written to disk.\n\n')
end
