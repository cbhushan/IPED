% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function estimate_3DSHORE_slice(DWI_file, bMatrices, user_options)
%  DWIfile - 4D nifti file name
%  bMatrices - 3 dimensional array of matrices (3x3xN)

disp('Estimating 3DSHORE ODFs in diffusion coordinates...')

% setting up the defaults options
opt = struct( ...
   'RadialOrder', 8, ...
   'shore_lambdaN', 1e-8, ...
   'shore_lambdaL', 1e-8, ...
   'shore_zeta', 700, ...
   'diffusion_time', 2.5*10^-3, ...
   'diffusion_coeff',0.7*10^-3,...
   'estimate_odf_3DSHORE', true, ...
   'SHORE_out_dir', fullfile(fileparts(DWI_file), '3DSHORE'), ...
   'bval_ratio_threshold', 45, ...
   'diffusion_modelling_mask', [], ...
   'Diffusion_coord_suffix', '.Diffusion.coord'... 
   );

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

% make the output directory
fname = fileBaseName(DWI_file);
SHORE_output_file_base = fullfile(opt.SHORE_out_dir, fname);
if opt.estimate_odf_3DSHORE && ~exist(opt.SHORE_out_dir, 'dir'), mkdir(opt.SHORE_out_dir); end

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
ind = 1:nDir; %find(~DEout.zero_bval_mask);
q = zeros(numel(ind),3);
del_t = opt.diffusion_time;
for i = 1:numel(ind)
    [V,D] = eig(bMatrices(:,:,ind(i)));
    [maxb,tmp] = max(diag(D));
    q(i,:) = V(:,tmp);
    qrad(i) = sqrt(maxb/(4*pi^2*del_t));
end

% compute mean b=0 image
if any(DEout.zero_bval_mask)
   b0mean = mean(dwi.img(:,:,:,DEout.zero_bval_mask), 4);
   zeroMask = zeroMask | (b0mean<=0);
else
   b0mean = ones(size(dwi.img,1), size(dwi.img,2), size(dwi.img,3));
end
%b0mean = b0mean(:);

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

outSize = [max(N), size(dwi.img,1), size(dwi.img,2), size(dwi.img,3)];

for nslice = 1:size(dwi.img,3),
   dwimages = double(reshape(dwi.img(:,:,nslice,:), [], nDir));
   dwimages = dwimages(:,ind);
   b0slice = vect(b0mean(:,:,nslice));
   dwimages = dwimages./b0slice(:,ones(1,size(dwimages,2))); % normalize by b=0 image
   dwimages(~isfinite(dwimages))=0;
   dwimages = permute(dwimages,[2,1]); % 1st dim is different diffusion weighting
   
   shore_coeff = shoreMatrixReg'*dwimages;
   
   % Convert to SH coefficients
   for num =  1:max(N),
      sh_3dshore(num,:,nslice) = sum(shore_coeff(N==num,:),1);
   end;
end;
clear dwi b0mean

sh_3dshore = permute(reshape(sh_3dshore,outSize), [2,3,4,1]);

shore_fid = fopen([SHORE_output_file_base '.SH.3DSHORE.odf'], 'w');

%GFA
% GFA = reshape(computeGFA_ODF_SHc(temp, 2*pi*Pf(L+1)), outSize(2:end));   
% temp3.img = GFA;
% temp3.img(zeroMask) = 0;
% fname = fullfile(fileparts(opt.FRT_out_dir), [fileBaseName(DWI_file) '.FRT_GFA.nii.gz']); % in parent folder of opt.FRT_out_dir
% save_untouch_nii_gz(temp3, fname, 16);

clear dwimages

% write output files
disp('Writing output 3DSHORE ODF files to disk...')
for k = 1:size(sh_3dshore,4) %Need to figure this out
    fname = sprintf('%s.SH.3DSHORE.%03d.nii.gz', SHORE_output_file_base, k);
    [~, name, ext] = fileparts(fname);
    fprintf(shore_fid, '%s\n', [name, ext]);
    temp3.img = sh_3dshore(:,:,:,k);
    temp3.img(zeroMask) = 0;
    save_untouch_nii_gz(temp3, fname, 16);
end
fclose(shore_fid);

fprintf('Estimated 3DSHORE ODFs written to disk.\n')

end
