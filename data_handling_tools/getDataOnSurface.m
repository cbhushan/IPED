% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [d] = getDataOnSurface(dfs, vol, method, str, ba_interp)
% Pulls out data at surface vertices/faces from corresponding 3D volume. dfs must be a surface
% saved by BrainSuite, corresponding to the 3D volume vol. 
% Usage:
%   d = getDataOnSurface(dfs, vol) % data on vertices
%   d = getDataOnSurface(dfs, vol, method) % defines interpolation method
%   d = getDataOnSurface(dfs, vol, method, 'face') % Gets data on faces
%   d = getDataOnSurface(dfs, vol, method, str, true) % Uses mex ba_interp3 (released under GPL)
%

if ~exist('ba_interp', 'var')
   ba_interp = false;
end

if ~exist('method', 'var')
   method = 'linear';
end

if ischar(vol)
   vol = load_untouch_nii_gz(vol, true);
elseif ~isfield(vol,'untouch') || vol.untouch ~= 1
   error('Please use ''load_untouch_nii.m'' to load the structure.');
end

if vol.hdr.hist.sform_code<=0
   error(['sform_code is not set in file %s. This file is not supported. You can try using '...
      'add_sform.m first.'], [fileBaseName(vol.fileprefix) '.nii']);
end

if ~isRAS_sform(vol)
   error('vol must be in RAS!')
end

if ndims(vol.img)<3
   error('vol must be atleast 3D volume.');
elseif ndims(vol.img)>4
   error('vol must be 3D or 4D volume!');
end

if ischar(dfs)
   dfs = readdfsGz(dfs);
elseif ~isfield(dfs,'vertices')
   error('dfs structure must have vertices field!');
end

res = vol.hdr.dime.pixdim(2:4);

if exist('str', 'var') && strcmpi(str, 'face')
   Fcntr = (dfs.vertices(dfs.faces(:,1),:) + dfs.vertices(dfs.faces(:,2),:) + dfs.vertices(dfs.faces(:,3),:))./3;
   Xtar = (Fcntr(:,1)./res(1))+1;
   Ytar = (Fcntr(:,2)./res(2))+1;
   Ztar = (Fcntr(:,3)./res(3))+1;
   
else % data on vertices
   Xtar = (dfs.vertices(:,1)./res(1))+1;
   Ytar = (dfs.vertices(:,2)./res(2))+1;
   Ztar = (dfs.vertices(:,3)./res(3))+1;
end


if ba_interp
   d = ba_interp3(double(vol.img), Ytar, Xtar, Ztar, method); % interpolates to border when outside
   d = squeeze(d);
   
else
   if ndims(vol.img)==3
      d = interpn(double(vol.img), Xtar, Ytar, Ztar, method, 0);      
   else
      nvol = size(vol.img,4);
      d = zeros(size(Xtar,1), nvol);
      for k = 1:nvol
         d(:,k) = interpn(double(vol.img(:,:,:,k)), Xtar, Ytar, Ztar, method, 0);
      end
   end
end

end

