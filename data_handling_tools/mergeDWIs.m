% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [fname, varargout] = mergeDWIs(out_file_base, varargin)
% Usage:
%   mergeDWIs(out_file_base, dwi_nifti_1, bval_1, bvec_1, dwi_nifti_2, bval_2, bvec_2, dwi_nifti_3, bval_3, bvec_3, ...)
%   mergeDWIs(out_file_base, dwi_nifti_1, bmat_1, dwi_nifti_2, bmat_2, dwi_nifti_3, bmat_3, ...)
%   mergeDWIs(out_file_base, dwi_nifti_1, bval_1, bvec_1, bmat_1, dwi_nifti_2, bval_2, bvec_2, bmat_2, dwi_nifti_3, bval_3, bvec_3, bmat_3, ...)

if nargin<3
   error('Incorrect usage. Number of inputs must be atleast 3. Check usage.')
end

opt_mode = parseinput(varargin{:});
fprintf('\nFiles per scan: %d', opt_mode);


% load first dataset
fprintf('\nReading dataset 1 ...\n')
dwi = load_untouch_nii(varargin{1});

if opt_mode==3 || opt_mode==4
   [bvec, bval] = readBvecBval(varargin{3}, varargin{2});
end

if opt_mode==2
   bMatrices = readBmat(varargin{2});
elseif opt_mode==4
   bMatrices = readBmat(varargin{4});
end

% concatenate data
c = 1;
for k = (1+opt_mode):opt_mode:length(varargin)
   c = c+1;
   fprintf('Reading dataset %d ...\n', c);
   temp = load_untouch_nii(varargin{k});
   dwi.img = cat(4, dwi.img, temp.img);
   
   if opt_mode==3 || opt_mode==4
      [temp1, temp2] = readBvecBval(varargin{k+2}, varargin{k+1});
      bvec = cat(1, bvec, temp1);
      bval = cat(1, bval, temp2);
   end
   
   if opt_mode==2
      temp = readBmat(varargin{k+1});
      bMatrices = cat(3, bMatrices, temp);
   elseif opt_mode==4
      temp = readBmat(varargin{k+3});
      bMatrices = cat(3, bMatrices, temp);
   end
end
clear temp*

% simple sanity check
if (opt_mode==3 || opt_mode==4) && (size(dwi.img, 4)~=length(bval) || size(bvec, 1)~=length(bval))
   error('Size mismatch in concatenated values. Make sure file names are in correct order!')
end

if (opt_mode==2 || opt_mode==4) && (size(dwi.img, 4)~=size(bMatrices, 3))
   error('Size mismatch in concatenated values. Make sure file names are in correct order!')
end

% save files
out_file_base = remove_extension(out_file_base);

fname = [out_file_base '.nii.gz'];
fprintf('Saving merged nifti file: %s\n', fname);
dwi.hdr.dime.dim(5) = size(dwi.img, 4);
save_untouch_nii_gz(dwi, fname);

varargout = {};
if (opt_mode==3 || opt_mode==4)
   bval_fname = [out_file_base '.bval'];
   bvec_fname = [out_file_base '.bvec'];
   fprintf('Saving merged bval/bvec file: \n%s\n%s\n', bval_fname, bvec_fname);
   writeBvecBvalFile(bvec, bvec_fname, bval, bval_fname);
   varargout{end+1} = bval_fname;
   varargout{end+1} = bvec_fname;
end

if (opt_mode==2 || opt_mode==4)
   bmat_fname = [out_file_base '.bmat'];
   fprintf('Saving merged bmat file: \n%s\n', bmat_fname);
   writeBmatFile(bMatrices, bmat_fname);
   varargout{end+1} = bmat_fname;
end

if nargout<1
   clear fname varargout
end
end

function n = parseinput(varargin)
if rem(nargin, 4)==0 && all(valid_img_filename(varargin(1:4:end))) % bval, bvec, bmat
   n = 4;
elseif rem(nargin, 3)==0 && all(valid_img_filename(varargin(1:3:end))) % bval, bvec
   n = 3;
elseif rem(nargin, 2)==0 && all(valid_img_filename(varargin(1:2:end))) % bmat
   n = 2;
else
   error('Something is wrong with inputs. Check number of inputs and filenames.')
end
end

function TF = valid_img_filename(fncell)
TF = false(length(fncell), 1);
for k = 1:length(fncell)
   TF(k) = check_endstr(fncell{k}, '.nii') | check_endstr(fncell{k}, '.nii.gz') | ...
      check_endstr(fncell{k}, '.img') | check_endstr(fncell{k}, '.hdr');
end
end

function TF = check_endstr(str, endstr)
loc = strfind(str, endstr);
len = length(str);
if (~isempty(loc)) && ((len-loc(end)) == (length(endstr)-1))
   TF = true;
else
   TF = false;
end
end

