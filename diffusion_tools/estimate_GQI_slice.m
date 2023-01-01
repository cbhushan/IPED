% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function estimate_GQI_slice(DWI_file, bMatrices, user_options)
%  DWIfile - 4D nifti file name
%  bMatrices - 3 dimensional array of matrices (3x3xN)

disp('Estimating GQI ODFs in diffusion coordinates...')

% setting up the defaults options
opt = struct( ...
   'HarmonicOrder', 8, ...
   'sigma_gqi', 1.25, ...
   'diffusion_time', 2.5*10^-3, ...
   'diffusion_coeff',0.7*10^-3,...
   'estimate_odf_GQI', true, ...
   'GQI_out_dir', fullfile(fileparts(DWI_file), 'GQI'), ...
   'bval_ratio_threshold', 45, ...
   'diffusion_modelling_mask', [], ...
   'Diffusion_coord_suffix', '.Diffusion.coord'... 
   );

% set options using user_options
if exist('user_options', 'var')
   if isfield(user_options, 'file_base_name')
      opt.GQI_out_dir = fullfile(fileparts(user_options.file_base_name), 'GQI');
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
GQI_output_file_base = fullfile(opt.GQI_out_dir, fname);
if opt.estimate_odf_GQI && ~exist(opt.GQI_out_dir, 'dir'), mkdir(opt.GQI_out_dir); end

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
    bval(i) = maxb;
end

% compute mean b=0 image
if any(DEout.zero_bval_mask)
   b0mean = mean(dwi.img(:,:,:,DEout.zero_bval_mask), 4);
   zeroMask = zeroMask | (b0mean<=0);
else
   b0mean = ones(size(dwi.img,1), size(dwi.img,2), size(dwi.img,3));
end
%b0mean = b0mean(:);

% GQI Computation
HarmonicOrder = opt.HarmonicOrder;
sigma_gqi = opt.sigma_gqi; % Take as input, for now hard coded
l_delta = sqrt(del_t*0.018); % 0.01506);  % diffusion length 0.01506 = sqrt(6*D), where D=2.5*10^-3 mm^2/sec % now hardcoded
odf_type = 1; % radial integral
lambda = 0.006;

% GQI  basis
basis = gqi_basis(sigma_gqi,l_delta,bval',[q(:,1),q(:,2),q(:,3)],[q(:,1),q(:,2),q(:,3)],odf_type,del_t);
[S,L] = sph_harm_basis([q(:,1),q(:,2),q(:,3)],HarmonicOrder,2);
sphericalHarmonicMatrixReg = (S'*S+lambda*diag(L.^2.*(L+1).^2))\(S');
gqiMatrixReg = sphericalHarmonicMatrixReg*basis;

outSize = [size(gqiMatrixReg,1), size(dwi.img,1), size(dwi.img,2), size(dwi.img,3)];

for nslice = 1:size(dwi.img,3),
   dwimages = double(reshape(dwi.img(:,:,nslice,:), [], nDir));
   dwimages = dwimages(:,ind);
   b0slice = vect(b0mean(:,:,nslice));
   dwimages = dwimages./b0slice(:,ones(1,size(dwimages,2))); % normalize by b=0 image
   dwimages(~isfinite(dwimages))=0;
   dwimages = permute(dwimages,[2,1]); % 1st dim is different diffusion weighting
   
   odf = basis'*dwimages;
   odf(odf<0) = 0;
   sh_gqi(:,:,nslice,:) = sphericalHarmonicMatrixReg*odf;
end;
clear dwi b0mean

sh_gqi = permute(reshape(sh_gqi,outSize), [2,3,4,1]);

gqi_fid = fopen([GQI_output_file_base '.SH.GQI.odf'], 'w');

%GFA
% GFA = reshape(computeGFA_ODF_SHc(temp, 2*pi*Pf(L+1)), outSize(2:end));   
% temp3.img = GFA;
% temp3.img(zeroMask) = 0;
% fname = fullfile(fileparts(opt.FRT_out_dir), [fileBaseName(DWI_file) '.FRT_GFA.nii.gz']); % in parent folder of opt.FRT_out_dir
% save_untouch_nii_gz(temp3, fname, 16);

clear dwimages

% write output files
disp('Writing output GQI ODF files to disk...')
for k = 1:size(sh_gqi,4) %Need to figure this out
    fname = sprintf('%s.SH.GQI.%03d.nii.gz', GQI_output_file_base, k);
    [~, name, ext] = fileparts(fname);
    fprintf(gqi_fid, '%s\n', [name, ext]);
    temp3.img = sh_gqi(:,:,:,k);
    temp3.img(zeroMask) = 0;
    save_untouch_nii_gz(temp3, fname, 16);
end
fclose(gqi_fid);

fprintf('Estimated GQI ODFs written to disk.\n')

end
