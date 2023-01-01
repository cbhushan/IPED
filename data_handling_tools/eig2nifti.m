% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [varargout] = eig2nifti(varargin)
%EIG2NIFTI convertes .eig.nii.gz files (as saved by BDP - BrainSuite Diffusion Pipeline) to standard
% NIfTI-1 format. 
%
% Usage:
%      eig2nifti <infname.eig.nii.gz> <outputfilebase>
%      eig2nifti(input_eig_file, output_base)
%
% input_eig_file is a character array of input filename (.eig.nii.gz) 
% output_base is character array of file-prefix for output filenames. 
% 
% Following example shows eig2nifti in action and explains the output files:
%
% >> eig2nifti('2467264c.dwi.RAS.correct.eig.nii.gz', '2467264c.tensor')
% Reading eig file...Done
% Saving files...
% Saved V1 file: 2467264c.tensor.V1.nii.gz
% Saved V2 file: 2467264c.tensor.V2.nii.gz
% Saved V3 file: 2467264c.tensor.V3.nii.gz
% Saved L1 file: 2467264c.tensor.L1.nii.gz
% Saved L2 file: 2467264c.tensor.L2.nii.gz
% Saved L3 file: 2467264c.tensor.L3.nii.gz
%
% L1 file contains the largest eigenvalues with corresponding eigenvector saved in V1 file. 
% Similarly, L2 file contains the second largest eigenvalues with corresponding eigenvector saved in
% V2 file followed by L3 file and V3 file. All eigenvalue files (L1,L2,L3) are 3D volume while 
% eigenvector files (V1,V2,V3) are 4D volume with last dimension containing the eigenvectors.
% 

if nargin ~= 2
   msg = '';
   if nargin ~= 0 
      msg = '\n\nIncorrect arguments. See Usage below.';
   end
   msg = {[msg '\n\neig2nifti convertes .eig.nii.gz files (as saved by BDP - BrainSuite '...
      'Diffusion Pipeline) to standard NIfTI-1 format. \n'], 'Usage:\n'...
      '\t eig2nifti <infname.eig.nii.gz> <outputfilebase>\n'...
      '\t eig2nifti(input_eig_file, output_base)\n\n',...
      'See usage in source file.'};
   error(sprintf(bdp_linewrap(msg)));
   
elseif ~ischar(varargin{1})
   error('The first argument should be a character array.')
   
elseif ~ischar(varargin{2})
   error('The second argument should be a character array.')
   
else
   input_eig_file = varargin{1};
   output_base =  remove_extension(varargin{2});
end

if exist(input_eig_file, 'file')~=2
   error('Specified eig file does not exist: %s', escape_filename(input_eig_file))
end

workdir = tempname;
mkdir(workdir);


% read .eig file
fprintf('\nReading eig file...')
eig = load_untouch_eig_gz(input_eig_file, workdir);
fprintf('Done')

% Eig structure 
% [v1.x v1.y v1.z v2.x v2.y v2.z v3.x v3.y v3.z l1 l2 l3]

fprintf('\nSaving files...')
temp4 = eig; 
temp4.img = [];
temp4.hdr.dime.dim(5) = 3;
temp4.hdr.hist.descrip = 'BrainSuite BDP - eig2nifti';

temp3 = temp4;
temp3.hdr.dime.dim(1) = 3;
temp3.hdr.dime.dim(5) = 1;

v1 = temp4;
v1.img = eig.img(:,:,:,1:3);
fname = save_untouch_nii_gz(v1, [output_base '.V1.nii.gz'], 16);
fprintf('\nSaved V1 file: %s', escape_filename(fname));

v2 = temp4;
v2.img = eig.img(:,:,:,4:6);
fname = save_untouch_nii_gz(v2, [output_base '.V2.nii.gz'], 16);
fprintf('\nSaved V2 file: %s', escape_filename(fname));

v3 = temp4;
v3.img = eig.img(:,:,:,7:9);
fname = save_untouch_nii_gz(v3, [output_base '.V3.nii.gz'], 16);
fprintf('\nSaved V3 file: %s', escape_filename(fname));

l1 = temp3;
l1.img = eig.img(:,:,:,10);
fname = save_untouch_nii_gz(l1, [output_base '.L1.nii.gz'], 16);
fprintf('\nSaved L1 file: %s', escape_filename(fname));

l2 = temp3;
l2.img = eig.img(:,:,:,11);
fname = save_untouch_nii_gz(l2, [output_base '.L2.nii.gz'], 16);
fprintf('\nSaved L2 file: %s', escape_filename(fname));

l3 = temp3;
l3.img = eig.img(:,:,:,12);
fname = save_untouch_nii_gz(l3, [output_base '.L3.nii.gz'], 16);
fprintf('\nSaved L3 file: %s\n\n', escape_filename(fname));

if nargout > 0
   varargout = {v1, v2, v3, l1, l2, l2};
end

end
