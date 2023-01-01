% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [data, X_vol, Y_vol, Z_vol, res] = get_grid_data(vol, method)
% Returns ndgrid type grid points and data (interpolated, if required) at those grid points. This
% can cause severe smoothing due to linear interpolation.
%

if ischar(vol)
   vol = load_untouch_nii_gz(vol);
elseif ~isfield(vol,'untouch') || vol.untouch ~= 1
   error('Please use ''load_untouch_nii.m'' to load the structure.');
end

if vol.hdr.hist.sform_code<=0 
   error('sform_code is not set in file %s. This file is not supported. You can try using add_sform.m first.', [file_base_name(vol.fileprefix) '.nii']);
end

% if vol.hdr.hist.qform_code>0 
%    fprintf('qform_code is set in %s. It is ignored and only sform matrix is used.\n', [file_base_name(vol.fileprefix) '.nii']);
% end


Tvol = zeros([4 4]);

Tvol(1,:)= vol.hdr.hist.srow_x(1:4);
Tvol(2,:)= vol.hdr.hist.srow_y(1:4);
Tvol(3,:)= vol.hdr.hist.srow_z(1:4);
Tvol(4,4) = 1;
t = Tvol(1:3, 1:3);

if isequal(diag(diag(t)), t)
   res = vol.hdr.dime.pixdim(2:4);
   data = vol.img;
   Tnew = Tvol;
else
   if ~exist('method', 'var') || strcmpi(method, 'linear')
      method = 1;
   elseif strcmpi(method, 'nearest')
      method = 2;
   else
      error('Unsupported method: %s', method);
   end
   
   res = min(vol.hdr.dime.pixdim(2:4))*[1 1 1];
   if ndims(vol.img) == 3
      [data Tnew] = affine(vol.img, Tvol, res, 0, 0, method);
      
   elseif ndims(vol.img) == 4
      [d Tnew] = affine(vol.img(:,:,:,1), Tvol, res, 0, 0, method);
      data = zeros([size(d), size(vol.img, 4)]);
      data(:,:,:,1) = d;
      for k = 2:size(vol.img, 4)
         data(:,:,:,k) = affine(vol.img(:,:,:,k), Tvol, res, 0, 0, method);
      end
   else
      error('vol has unsupported dimensions.');
   end
end

clear vol

vol_size = size(data);
[X_vol, Y_vol, Z_vol] = ndgrid(0:vol_size(1)-1, 0:vol_size(2)-1,  0:vol_size(3)-1);
c(1,:) = X_vol(:); clear X_vol;
c(2,:) = Y_vol(:); clear Y_vol;
c(3,:) = Z_vol(:);
c(4,:) = 1;
c = Tnew*c;
X_vol = reshape(c(1,:), size(Z_vol));
Y_vol = reshape(c(2,:), size(Z_vol));
Z_vol = reshape(c(3,:), size(Z_vol));
clear c

end
