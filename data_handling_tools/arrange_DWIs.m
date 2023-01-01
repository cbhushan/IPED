% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [dwi, b0s, mean_b0, varargout] = arrange_DWIs(DWI_file, varargin)
% Rearranges DWIs (& bmat/bval/bvec) files so that all b=0 images (& corresponding
% bval/bvec/bmat) come before diffusion weighted images. Writes out re-arranged files to disk.
%
% Usage: 
%    [dwi_a, b0s, mean_b0, bmat_a] = arrange_DWIs(DWI_file, bmat_file, output_file_base)
%    [dwi_a, b0s, mean_b0, bval_a, bvec_a] = arrange_DWIs(DWI_file, bval_file, bvec_file, output_file_base)
%    [dwi_a, b0s, mean_b0, bval_a, bvec_a, bmat_a] = arrange_DWIs(DWI_file, bval_file, bvec_file, bmat_file, output_file_base)
%

workDir = tempname;
mkdir(workDir);

% parse input
if length(varargin) < 2
   error('Too few input.');
   
elseif length(varargin) == 2   % bmat
   bmat_file = varargin{1};
   
elseif length(varargin) == 3   % bval, bvec
   bval_file = varargin{1};
   bvec_file = varargin{2};

elseif length(varargin) == 4  % bval, bvec, bmat
   bval_file = varargin{1};
   bvec_file = varargin{2};
   bmat_file = varargin{3};
end
output_file_base = varargin{end};


% load bval, bvec, bmat
if exist('bval_file', 'var')
   [bvec, bval] = readBvecBval(bvec_file, bval_file);
   bval2use = bval;
end

if exist('bmat_file', 'var')
   bMatrices = readBmat(bmat_file);
   nDir = size(bMatrices, 3);
   bval2use = zeros(1,nDir);
   for iDir = 1:nDir
      bval2use(iDir) = trace(bMatrices(:,:,iDir));
   end
end

% sort by b-value
[temp, Ix] = sort(bval2use, 'ascend');
varargout = {};

% re-arrange bvals
if exist('bval', 'var')
   bval = bval(Ix);
   bvec = bvec(Ix,:);
   writeBvecBvalFile(bvec, [remove_extension(output_file_base) '.bvec'], ...
      bval, [remove_extension(output_file_base) '.bval'])
   varargout{end+1} = bval;
   varargout{end+1} = bvec;
   DE_out = checkDiffusionEncodingScheme(bvec, bval);
end

% re-arrange bmatrices
if exist('bMatrices', 'var')
   bMatrices = bMatrices(:,:,Ix);
   writeBmatFile(bMatrices, [remove_extension(output_file_base) '.bmat'])
   varargout{end+1} = bMatrices;
   DE_out = checkDiffusionEncodingScheme(bMatrices);
end

% re-arrange dwi
dwi = load_untouch_nii_gz(DWI_file, workDir);
dwi.img = dwi.img(:,:,:,Ix);
save_untouch_nii_gz(dwi, remove_extension(output_file_base), workDir);

num_b0s = sum(DE_out.zero_bval_mask);
b0s = dwi;
b0s.img = dwi.img(:,:,:,DE_out.zero_bval_mask);
b0s.hdr.dime.dim(5) = num_b0s;
save_untouch_nii_gz(b0s, [remove_extension(output_file_base) '.b0s'], workDir);

mean_b0 = b0s;
mean_b0.img = mean(double(dwi.img(:,:,:,DE_out.zero_bval_mask)), 4);
mean_b0.hdr.dime.dim(5) = 1;
mean_b0.hdr.dime.dim(1) = 3;
save_untouch_nii_gz(mean_b0, [remove_extension(output_file_base) '.mean_b0'], workDir);

rmdir(workDir, 's')
end
