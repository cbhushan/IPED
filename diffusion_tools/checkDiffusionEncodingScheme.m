% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function out = checkDiffusionEncodingScheme(varargin)
% Checks few conditions diffusion encoding scheme. See description of out below. 
%
% Usage:
%    out = checkDiffusionEncodingScheme(bvec, bval)
%    out = checkDiffusionEncodingScheme(bmat)
%    out = checkDiffusionEncodingScheme(bvec, bval, bval_thresh)
%    out = checkDiffusionEncodingScheme(bmat, bval_thresh)
%
% out.bCond - condition number for un-weighted tensor fit
% out.illposed_tensor_fit - boolean, true when illposed 
% out.zero_bval_mask - estimated b=0 mask (may not be pefect)
% out.bval - estimated b-values
%

bval_thresh = 45; % ratio of 900 to 20
bval_thresh_nz = 1.5; % 1500/1000

if nargin==1 || (nargin==2 && numel(varargin{2})==1 && isnumeric(varargin{2}))
   if ischar(varargin{1})
      bMatrices = readBmat(varargin{1});
   else
      bMatrices = varargin{1};
   end   
   [bval, bvec] = bmat2BvecBval(bMatrices);

   if nargin==2
      bval_thresh = varargin{2};
   end   
   
elseif nargin==3 || nargin==2
   if ischar(varargin{1}) && ischar(varargin{2})
      [bvec, bval] = readBvecBval(varargin{1}, varargin{2});
   else
      bvec = varargin{1};
      bval = varargin{2};
      
      if size(bvec,1)~=length(bval)
         err_msg = {'Number of values in bvec and bval do not match.', ...
            ['Please make sure that they have equal number of entries corresponding '...
            'to diffusion encoding directions.']};
         error('BDP:InconsistentDiffusionParameters', bdp_linewrap(err_msg))
      end
   end
   
   % check bvec for unit Norm
   bvec_norm = sqrt(sum(bvec.^2, 2));
   bvec_norm(bvec_norm<0.01) = []; % remove norm corresponding to b=0
   if any(abs(bvec_norm-1)>0.01)
      err_msg = {['Few diffusion encoding directions in bvec file are not unit norm. ', ...
         'All diffusion gradient vectors must be unit norm.']};
      error('BDP:InvalidFile', bdp_linewrap(err_msg));
   end
   
   vec = bvec';
   bMatrices = zeros(3,3,size(vec,2));
   for iDir = 1:size(vec,2)
      bMatrices(:,:,iDir) = bval(iDir)*( vec(:,iDir)*(vec(:,iDir)') );
   end
   
   if nargin==3
      bval_thresh = varargin{3};
   end   
else
   error('Number of inputs should be exactly 1, 2 or 3')
end

% add parameters to output
out.bval = bval;
out.bvec = bvec;
out.bmat = bMatrices;


% Check for ill-posed tensor fit
nDir = size(bMatrices, 3);
bMatrix = zeros(nDir,7);
for bIndex = 1:nDir
    bMatrix(bIndex,:) = -[bMatrices(1,1,bIndex) bMatrices(2,2,bIndex) bMatrices(3,3,bIndex) ...
                          2*bMatrices(2,1,bIndex) 2*bMatrices(3,1,bIndex) 2*bMatrices(3,2,bIndex) -1];
end
out.bCond = cond(bMatrix, 2);
out.illposed_tensor_fit = out.bCond>1e6;

% Try to find b=0. Very simple 
out.zero_bval_mask = false(size(bval));
if min(bval)==0
   out.zero_bval_mask = (bval==0);
else   
   b_low = min(bval);
   b_high = max(bval);
   if (b_high/b_low)>=bval_thresh 
      out.zero_bval_mask = (bval==b_low);      
   end
end

% Find the range of the remaining b-values
bval_nzero = bval(out.zero_bval_mask==0);
b_low = min(bval_nzero);
b_high = max(bval_nzero);
if (b_high/b_low)>=bval_thresh_nz 
   out.single_shell =0;
else
   out.single_shell =1;
end;

end

function [bval, bvec] = bmat2BvecBval(bMatrices)
bval = zeros(size(bMatrices, 3), 1);
bvec = zeros(size(bMatrices, 3), 3);
for iDir = 1:size(bMatrices, 3)
   [V, D] = sorted_eig(bMatrices(:,:,iDir));
   bval(iDir) = D(1);
   bvec(iDir,:) = transpose(V(:,1));
end
end

