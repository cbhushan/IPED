% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function estimate_tensors(DWIfile, bMatrices, user_options)
%  Estimates diffusion tensors from DWIs & saves files in nifti format. 
%  Weighted Log Linear fit. No eddy current correction or motion correction.
%
%  Assumes bMatrices to be specified in voxel-coordinates of DWIfile.
%  Please input correctly rotated bMatrices.  
%
%  DWIfile - 4D nifti file name
%  bMatrices - 3 dimensional array of matrices (3x3xN)

fprintf('\nEstimating diffusion tensors in diffusion coordinates...')

opt = struct( ...
   'tensor_out_dir', fileparts(DWIfile), ...
   'bval_ratio_threshold', 45, ...
   'diffusion_modelling_mask', [], ...
   'Diffusion_coord_suffix', '.Diffusion.coord'...
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
output_file_base = fullfile(opt.tensor_out_dir, fileBaseName(DWIfile));

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

dwi = load_untouch_nii_gz(DWIfile, true);
nDir = size(dwi.img, 4);

if nDir~=size(bMatrices, 3)
   error('BDP:InconsistentDiffusionParameters', ['Number of diffusion image should be exactly same' ...
      'as number of bMatrices/bvec (in .bmat/.bvec file, if any).']);
end

% load zero mask
if isempty(opt.diffusion_modelling_mask)
   sz = size(dwi.img);
   zeroMask = false(sz(1:3));
elseif exist(opt.diffusion_modelling_mask, 'file')~=2
   error(['DWI mask not found: ' escape_filename(opt.diffusion_modelling_mask)])
else
   b0s = dwi;
   b0s.img = dwi.img(:,:,:,DEout.zero_bval_mask);
   b0s.hdr.dime.dim(5) = sum(DEout.zero_bval_mask);
   fname = save_untouch_nii_wrapper(b0s, [tempname() '.nii']);
   temp = fixBSheader(fname, opt.diffusion_modelling_mask);
   zeroMask = temp.img<=0;
   clear b0s temp
end

% remove voxels with zeros in b=0 images
zeroMask = zeroMask | ~all(dwi.img(:,:,:, DEout.zero_bval_mask), 4);

% Removes voxels in diffusion-coord mask which are part of zeroMask - SHADY CODE, Maybe FIX IT
if exist('user_options', 'var') && isfield(user_options, 'file_base_name')...
      && isfield(user_options, 't1_mask_file') && ~isempty(user_options.t1_mask_file)
   
   [~, nm, ext] = fileparts(user_options.t1_mask_file);
   msk_fname = fullfile(fileparts(user_options.file_base_name), ...
      suffix_filename([nm ext], user_options.Diffusion_coord_suffix));
   
   if exist(msk_fname, 'file')
      temp = load_untouch_nii_gz(msk_fname);
      temp.img(zeroMask) = 0;
      save_untouch_nii_gz(temp, msk_fname);
   end   
end

% fit diffusion tensors
temp = double(dwi.img(~zeroMask(:,:,:,ones(1,nDir))));
temp = permute(reshape(temp, [], nDir), [2,1]); % 1st dim is different diffusion weighting
[lmbd, eigVec] = dtiLogLinearW(temp, bMatrices);
mADC_vec = compute_mADC(temp, DEout);
clear temp

lambdas  = zeros([size(zeroMask) 3]);
lambdas(~zeroMask(:,:,:,[1 1 1])) = vect(permute(lmbd,[2 1]));

eigenVectors = zeros([size(zeroMask) 3 3]);
eigenVectors(~zeroMask(:,:,:,[1 1 1],[1 1 1])) = vect(permute(eigVec,[3 1 2]));

lambdas = abs(permute(lambdas, [4 1 2 3]));
eigenVectors = permute(eigenVectors, [4 5 1 2 3]);

mADC = zeros(size(zeroMask));
mADC(~zeroMask) = mADC_vec(:);
clear mADC_vec

fprintf('\nWriting tensor files to disk...')

% template 3D volume
temp3 = dwi;
temp3.img = 0;
temp3.hdr.dime.dim(1) = 3;
temp3.hdr.dime.dim(5) = 1;

% template 4D volume
temp4 = dwi;
temp4.img = 0;
temp4.hdr.dime.dim(5) = 3;

v1 = temp4;
v1.hdr.dime.dim(5) = 3;
v1.img = squeeze(eigenVectors(:,1,:,:,:));
v1.img = permute(v1.img, [2 3 4 1]);
v1.img(repmat(zeroMask, [1,1,1,3])) = 0;

v2 = temp4;
v2.img = squeeze(eigenVectors(:,2,:,:,:));
v2.img = permute(v2.img, [2 3 4 1]);
v2.img(repmat(zeroMask, [1,1,1,3])) = 0;

v3 = temp4;
v3.img = squeeze(eigenVectors(:,3,:,:,:));
v3.img = permute(v3.img, [2 3 4 1]);
v3.img(repmat(zeroMask, [1,1,1,3])) = 0;

FA = temp3;
den = 2*( lambdas(1,:,:,:).^2 + lambdas(2,:,:,:).^2 + lambdas(3,:,:,:).^2 );
num = (lambdas(1,:,:,:)-lambdas(2,:,:,:)).^2 ...
     +(lambdas(2,:,:,:)-lambdas(3,:,:,:)).^2 ...
     +(lambdas(3,:,:,:)-lambdas(1,:,:,:)).^2;
FA.img = squeeze(sqrt(num./den));
zmask = (squeeze(den)==0);
FA.img(zeroMask | zmask) = 0;
save_untouch_nii_gz(FA,[output_file_base '.FA.nii.gz'], 16);

colorFA = temp3;
colorFA.img = zeros([size(FA.img) 3]);
colorFA.img(:,:,:,1) = abs(squeeze(eigenVectors(1,1,:,:,:)).*FA.img)*255;
colorFA.img(:,:,:,2) = abs(squeeze(eigenVectors(2,1,:,:,:)).*FA.img)*255;
colorFA.img(:,:,:,3) = abs(squeeze(eigenVectors(3,1,:,:,:)).*FA.img)*255;
colorFA.img(colorFA.img>255) = 0;
save_untouch_nii_gz(colorFA,[output_file_base '.FA.color.nii.gz'], 128);
clear colorFA FA eigenVectors


L1 = temp3;
L1.img = squeeze(lambdas(1,:,:,:));
L1.img(zeroMask) = 0;

L2 = temp3;
L2.img = squeeze(lambdas(2,:,:,:));
L2.img(zeroMask) = 0;
save_untouch_nii_gz(L2, [output_file_base '.L2.nii.gz'], 16);

L3 = temp3;
L3.img = squeeze(lambdas(3,:,:,:));
L3.img(zeroMask) = 0;
save_untouch_nii_gz(L3, [output_file_base '.L3.nii.gz'], 16);

axial = temp3;
axial.img = squeeze(lambdas(1,:,:,:));
save_untouch_nii_gz(axial, [output_file_base '.axial.nii.gz'], 16);
clear axial

radial = temp3;
radial.img = squeeze(mean(lambdas(2:3,:,:,:), 1));
save_untouch_nii_gz(radial, [output_file_base '.radial.nii.gz'], 16);
clear radial

MD = temp3;
MD.img = squeeze(mean(lambdas,1));
save_untouch_nii_gz(MD, [output_file_base '.MD.nii.gz'], 16);
clear MD lambdas

temp3.img = mADC;
save_untouch_nii_gz(temp3, [output_file_base '.mADC.nii.gz'], 16);
clear temp3 mADC

generate_eig_file(L1, L2, L3, v1, v2, v3, [output_file_base '.eig.nii.gz'] );

fprintf('\nEstimated diffusion tensors parameters written to disk.\n')
end
