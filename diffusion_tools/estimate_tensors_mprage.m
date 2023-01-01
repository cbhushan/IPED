% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function estimate_tensors_mprage(data_file, target_file, affinematrixfile, bMatrices, user_options)
%  Estimates Diffusion Tensors from Diffusion images & saves files in nifti format. Saves the output in
% BrainSuite friendly format. 
%
%  No eddy current or motion correction.
%
%  data_file - 4D nifti diffusion file name
%  target_file - bfc file name (this will be also used as mask)
%  affinematrixfile - file name of .mat file saved after rigid registration of diffusion and
%                     bfc file
%  bMatrices - .bmat file name or array of matrices (3x3xN)
%

fprintf('\nEstimating diffusion tensors in T1 coordinate...')

opt = struct( ...
   'tensor_out_dir', fileparts(data_file), ...
   'mprage_coord_suffix', '.T1_coord', ...
   'bval_ratio_threshold', 45, ...
   't1_mask_file', target_file ...
   );

% set options using user_options
if exist('user_options', 'var')
   if isfield(user_options, 'file_base_name')
      opt.tensor_out_dir = fileparts(user_options.file_base_name);
   end
   
   fnames = fieldnames(opt);
   for iname = 1:length(fnames)
      if isfield(user_options, fnames{iname})
         opt = setfield(opt, fnames{iname}, getfield(user_options, fnames{iname}));
      end
   end
end
output_file_base = fullfile(opt.tensor_out_dir, fileBaseName(data_file));

% Apply the rigid transform to header of dwi file
load(affinematrixfile);
data_file_transformed = [tempname '.nii.gz'];
[~,Tnew] = affine_transform_nii(data_file, M_world, origin, data_file_transformed);
clear Tnew;

% load bmatrices
if ischar(bMatrices)
   bMatrices = readBmat(bMatrices);
end

DEout = checkDiffusionEncodingScheme(bMatrices, opt.bval_ratio_threshold);
if DEout.illposed_tensor_fit
   warn_msg = ['Based on the the input diffusion encoding parameters (bval, bvec, bmat), the '...
      'diffusion tensor (DT) fit could be ill-conditioned. This could be result of missing acquisition '...
      'of the diffusion image corresponding to b-value of 0 (zero). OR, it could be an effect of insufficient '...
      'distinct diffusion encoding directions. However, BDP will try to estimate diffusion tensors. Please '...
      'check the outputs and be careful with the analysis of the estimated diffusion tensors.'];
   bdpPrintWarning('Ill-conditioned DT fit', warn_msg);
end
if ~DEout.single_shell
     msg = 'The input diffusion data seems to have a non-single-shell q-space sampling pattern. DTI is a single-shell sampling based method and can give unexpected tensor results when used with differently sampled diffusion MRI data. Please be careful with interpretation and analysis of the outputs. This warning can be ignored if the input data was indeed sampled with single-shell acquisition.';
     bdpPrintWarning('Inappropriate diffusion model selected:', msg);
end;
 
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

% make mask
vd = load_untouch_nii_gz(opt.t1_mask_file);
targetMask = vd.img>0;
vd = load_untouch_nii_gz(target_file,true);
vd.img = [];

% template 3D volume
temp3 = vd;
temp3.img = 0;
temp3.hdr.dime.dim(1) = 3;
temp3.hdr.dime.dim(5) = 1;

% template 4D volume
temp4 = vd;
temp4.img = 0;
temp4.hdr.dime.dim(1) = 4;
temp4.hdr.dime.dim(5) = 3;


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

lambdas = zeros([3 target_vol_size]);
eigenVectors = zeros([3 3 target_vol_size]);
mADC = zeros(target_vol_size);
cpb = ConsoleProgressBar(); % Set progress bar parameters
cpb.setMinimum(1);
cpb.setMaximum(target_vol_size(3));
cpb.start();
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
      zeroMask = ~all(dwimages(:,DEout.zero_bval_mask), 2); % voxels with any zero b=0 image
      
      dwimages = permute(dwimages,[2,1]);
      zeroMask = zeroMask';
      
      % fit diffusion tensors
      [lbd, eigVec] = dtiLogLinearW(dwimages, bMatrices);
      lbd(zeroMask([1 1 1], :)) = 0;
      mADC_vec = compute_mADC(dwimages, DEout);
      mADC_vec(zeroMask) = 0;
      
      msk4D = permute(targetMask(:,:,slice, [1 1 1]), [4 1 2 3]);
      outimg = zeros([3 target_vol_size(1:2) 1]);
      outimg(msk4D) = lbd(:);
      lambdas(:,:,:,slice) = outimg;
      
      msk4D = permute(targetMask(:,:,slice, [1 1 1], [1 1 1]), [4 5 1 2 3]);
      outimg = zeros([3 3 target_vol_size(1:2) 1]);
      outimg(msk4D) = eigVec(:);
      eigenVectors(:,:,:,:,slice) = outimg;
      
      outimg = zeros([target_vol_size(1:2) 1]);
      outimg(msk) = mADC_vec(:);
      mADC(:,:,slice) = outimg;
   end
   
   text = sprintf('%d/%d slices done', slice, target_vol_size(3));
   cpb.setValue(slice); cpb.setText(text);
end
cpb.stop();
lambdas = abs(lambdas);
fprintf('\nWriting tensor files to disk...')


v1 = temp4;
v1.hdr.dime.dim(5) = 3;
v1.img = squeeze(eigenVectors(:,1,:,:,:));
v1.img = permute(v1.img, [2 3 4 1]);

v2 = temp4;
v2.img = squeeze(eigenVectors(:,2,:,:,:));
v2.img = permute(v2.img, [2 3 4 1]);

v3 = temp4;
v3.img = squeeze(eigenVectors(:,3,:,:,:));
v3.img = permute(v3.img, [2 3 4 1]);

L1 = temp3;
L1.img = squeeze(lambdas(1,:,:,:));

L2 = temp3;
L2.img = squeeze(lambdas(2,:,:,:));
save_untouch_nii_gz(L2,[output_file_base '.L2' opt.mprage_coord_suffix '.nii.gz'], 16);

L3 = temp3;
L3.img = squeeze(lambdas(3,:,:,:));
save_untouch_nii_gz(L3,[output_file_base '.L3' opt.mprage_coord_suffix '.nii.gz'], 16);

temp = temp3;
temp.img = mADC;
save_untouch_nii_gz(temp ,[output_file_base '.mADC' opt.mprage_coord_suffix '.nii.gz'], 16);
clear temp mADC

generate_eig_file(L1, L2, L3, v1, v2, v3, [output_file_base opt.mprage_coord_suffix '.eig.nii.gz']);
clear L1 L2 L3 v1 v2 v3;

axial = temp3;
axial.img = squeeze(lambdas(1,:,:,:));
save_untouch_nii_gz(axial,[output_file_base '.axial' opt.mprage_coord_suffix '.nii.gz'], 16);
clear axial

radial = temp3;
radial.img = squeeze(mean(lambdas(2:3,:,:,:),1));
save_untouch_nii_gz(radial,[output_file_base '.radial' opt.mprage_coord_suffix '.nii.gz'], 16);
clear radial

MD = temp3;
MD.img = squeeze(mean(lambdas,1));
save_untouch_nii_gz(MD,[output_file_base '.MD' opt.mprage_coord_suffix '.nii.gz'], 16);
clear MD

FA = temp3;
den = 2*( lambdas(1,:,:,:).^2 + lambdas(2,:,:,:).^2 + lambdas(3,:,:,:).^2 );
zmsk = (den==0);
num = (lambdas(1,:,:,:)-lambdas(2,:,:,:)).^2 ...
     +(lambdas(2,:,:,:)-lambdas(3,:,:,:)).^2 ...
     +(lambdas(3,:,:,:)-lambdas(1,:,:,:)).^2;
FA.img = squeeze(sqrt(num./den));
FA.img(zmsk) = 0;
save_untouch_nii_gz(FA,[output_file_base '.FA' opt.mprage_coord_suffix '.nii.gz'], 16);


colorFA = temp3;
colorFA.img = zeros([size(FA.img) 3]);
colorFA.img(:,:,:,1) = abs(squeeze(eigenVectors(1,1,:,:,:)).*FA.img)*255;
colorFA.img(:,:,:,2) = abs(squeeze(eigenVectors(2,1,:,:,:)).*FA.img)*255;
colorFA.img(:,:,:,3) = abs(squeeze(eigenVectors(3,1,:,:,:)).*FA.img)*255;
colorFA.img(colorFA.img>255) = 0;
save_untouch_nii_gz(colorFA,[output_file_base '.FA.color' opt.mprage_coord_suffix '.nii.gz'], 128);
clear colorFA FA;


if exist(data_file_transformed,'file')
   delete(data_file_transformed);
end

fprintf('\nEstimated diffusion tensors parameters written to disk.\n')
end
