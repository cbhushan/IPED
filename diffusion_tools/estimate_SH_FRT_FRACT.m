% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function estimate_SH_FRT_FRACT(DWI_file, bMatrices, user_options)
%  DWIfile - 4D nifti file name
%  bMatrices - 3 dimensional array of matrices (3x3xN)

disp('Estimating ODFs in diffusion coordinates...')

% setting up the defaults options
opt = struct( ...
   'HarmonicOrder', 8, ...
   'xi', 0.34, ...
   'odf_lambda', 0.006, ...
   'FRT_out_dir', fullfile(fileparts(DWI_file), 'FRT'), ...
   'FRACT_out_dir', fullfile(fileparts(DWI_file), 'FRACT'), ...
   'estimate_odf_FRT', false, ...
   'estimate_odf_FRACT', true, ...
   'bval_ratio_threshold', 45, ...
   'diffusion_modelling_mask', [], ...
   'Diffusion_coord_suffix', '.Diffusion.coord'... 
   );

% set options using user_options
if exist('user_options', 'var')
   if isfield(user_options, 'file_base_name')
      opt.FRT_out_dir = fullfile(fileparts(user_options.file_base_name), 'FRT');
      opt.FRACT_out_dir = fullfile(fileparts(user_options.file_base_name), 'FRACT');
   end
   
   fnames = fieldnames(opt);
   for iname = 1:length(fnames)
      if isfield(user_options, fnames{iname})
         opt = setfield(opt, fnames{iname}, getfield(user_options, fnames{iname}));
      end
   end
end

fname = fileBaseName(DWI_file);
FRT_output_file_base = fullfile(opt.FRT_out_dir, fname);
FRACT_output_file_base = fullfile(opt.FRACT_out_dir, fname);
if opt.estimate_odf_FRT && ~exist(opt.FRT_out_dir, 'dir'), mkdir(opt.FRT_out_dir); end
if opt.estimate_odf_FRACT && ~exist(opt.FRACT_out_dir, 'dir'), mkdir(opt.FRACT_out_dir); end

% load data
dwi = load_untouch_nii_gz(DWI_file, true);
nDir = size(dwi.img, 4);

% load bmatrices
if ischar(bMatrices)
   bMatrices = readBmat(bMatrices);
end
DEout = checkDiffusionEncodingScheme(bMatrices, opt.bval_ratio_threshold);

if nDir~=size(bMatrices, 3)
   error('BDP:InconsistentDiffusionParameters', ['Number of diffusion image should be exactly same' ...
      'as number of bMatrices/bvec (in .bmat/.bvec file, if any).']);
end

if ~DEout.single_shell
     msg = 'The input diffusion data seems to have a non-single-shell q-space sampling pattern. FRT and FRACT are single-shell based ODF estimation methods and can give unexpected ODF results when used with differently sampled diffusion MRI data. Please be careful with interpretation and analysis of the outputs. This warning can be ignored if the input data was indeed sampled with single-shell acquisition.';
     bdpPrintWarning('Inappropriate diffusion model selected:', msg);
 end;
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

% template 3D volume
temp3 = dwi;
temp3.img = 0;
temp3.hdr.dime.dim(1) = 3;
temp3.hdr.dime.dim(5) = 1;

% Find diffusion encoding direction (same as largest eigen vector of bMatrices)
ind = find(~DEout.zero_bval_mask);
q = zeros(numel(ind),3);
for i = 1:numel(ind)
    [V,D] = eig(bMatrices(:,:,ind(i)));
    [~,tmp] = max(diag(D));
    q(i,:) = V(:,tmp);
end

% compute mean b=0 image
if any(DEout.zero_bval_mask)
   b0mean = mean(dwi.img(:,:,:,DEout.zero_bval_mask), 4);
   zeroMask = zeroMask | (b0mean<=0);
else
   b0mean = ones(size(dwi.img,1), size(dwi.img,2), size(dwi.img,3));
end
b0mean = b0mean(:);

% FRT and FRACT Computation
HarmonicOrder = opt.HarmonicOrder;
xi = opt.xi;
lambda = opt.odf_lambda;

[S,L] = sph_harm_basis([q(:,1),q(:,2),q(:,3)],HarmonicOrder,2);

for order = 0:2:HarmonicOrder
    lpp = legendre(order,xi);
    lpo = legendre(order,0);
    lpn = legendre(order,-xi);
    Pg(order+1) = 2*lpo(1)-lpp(1)-lpn(1);
    Pf(order+1) = 2*lpo(1);
end

sphericalHarmonicMatrixReg = (S'*S+lambda*diag(L.^2.*(L+1).^2))\(S');

dwimages = double(reshape(dwi.img, [], nDir));
dwimages = dwimages(:,ind);
dwimages = dwimages./b0mean(:,ones(1,size(dwimages,2))); % normalize by b=0 image
dwimages(~isfinite(dwimages))=0;
dwimages = permute(dwimages,[2,1]); % 1st dim is different diffusion weighting
outSize = [size(S,2), size(dwi.img,1), size(dwi.img,2), size(dwi.img,3)];
clear dwi b0mean

if opt.estimate_odf_FRACT
   fractMatrixReg = diag(2*pi*Pg(L+1))*sphericalHarmonicMatrixReg;
   sh_fract = permute(reshape(fractMatrixReg*dwimages, outSize), [2,3,4,1]);
   fract_fid = fopen([FRACT_output_file_base '.SH.FRACT.odf'], 'w');
end

if opt.estimate_odf_FRT
   frtMatrixReg = diag(2*pi*Pf(L+1))*sphericalHarmonicMatrixReg;
   temp = frtMatrixReg*dwimages;
   sh_frt = permute(reshape(temp, outSize), [2,3,4,1]);
   frt_fid = fopen([FRT_output_file_base '.SH.FRT.odf'], 'w');
   
   GFA = reshape(computeGFA_ODF_SHc(temp, 2*pi*Pf(L+1)), outSize(2:end));   
   temp3.img = GFA;
   temp3.img(zeroMask) = 0;
   fname = fullfile(fileparts(opt.FRT_out_dir), [fileBaseName(DWI_file) '.FRT_GFA.nii.gz']); % in parent folder of opt.FRT_out_dir
   save_untouch_nii_gz(temp3, fname, 16);
end
clear dwimages


% write output files
disp('Writing output ODF files to disk...')
for k = 1:size(S,2)
   if opt.estimate_odf_FRACT
      fname = sprintf('%s.SH.FRACT.%03d.nii.gz', FRACT_output_file_base, k);
      [~, name, ext] = fileparts(fname);
      fprintf(fract_fid, '%s\n', [name, ext]);
      temp3.img = sh_fract(:,:,:,k);
      temp3.img(zeroMask) = 0;
      save_untouch_nii_gz(temp3, fname, 16);
   end
   
   if opt.estimate_odf_FRT
      fname = sprintf('%s.SH.FRT.%03d.nii.gz', FRT_output_file_base, k);
      [~, name, ext] = fileparts(fname);
      fprintf(frt_fid, '%s\n', [name, ext]);
      temp3.img = sh_frt(:,:,:,k);
      temp3.img(zeroMask) = 0;
      save_untouch_nii_gz(temp3, fname, 16);
   end
end

if opt.estimate_odf_FRACT, fclose(fract_fid); end
if opt.estimate_odf_FRT, fclose(frt_fid); end

fprintf('Estimated ODFs written to disk.\n')

end

