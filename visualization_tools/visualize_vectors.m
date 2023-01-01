% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function visualize_vectors(nifti_file, png_file, bg_file)

close all
v1 = load_untouch_nii_gz(nifti_file);
if exist('bg_file', 'var')
   bg = load_untouch_nii_gz(bg_file);
else
   bg.img = zeros(size(v1.img));
end

figSize = 1050;

% sform
sformT = zeros([3 3]);
sformT(1,:)= v1.hdr.hist.srow_x(1:3);
sformT(2,:)= v1.hdr.hist.srow_y(1:3);
sformT(3,:)= v1.hdr.hist.srow_z(1:3);

if ~isequal(sformT, diag(diag(abs(sformT))))
   sformT
   warning('sformT is not +ve definite diagonal matrix. It will be ignored.')
end

v1size = size(v1.img);
if length(v1size)<4
   error('v1 should be 4D matrix');
end

res = v1.hdr.dime.pixdim(2:4);

% normalize
v1.img(~isfinite(v1.img)) = 0;
vec_len = sqrt(sum(v1.img.^2, 4));
v1.img = v1.img ./ vec_len(1:end, 1:end, 1:end, [1 1 1]);
v1.img(~isfinite(v1.img)) = 0;

[xgrid ygrid] = ndgrid(res(1):res(1):v1size(1)*res(1), res(2):res(2):v1size(2)*res(2));

hfig = figure('units','pixels', 'position', [20, 20, figSize, figSize]); %'color', [0 0 0]);
haxes = axes('units','normalized', 'position', [0, 0, 1, 1], 'color', [0 0 0]);
hold on

for slice_no = round(v1size(3)/2-5):round(v1size(3)/2+5) %30:401:v1size(3)
   bg_img = squeeze(bg.img(:,:,slice_no))';
   imagesc([res(1) v1size(1)*res(1)], [res(2) v1size(2)*res(2)], bg_img, [0 1]);
   colormap(gray)
   set(gca,'ydir','normal');
   hold on;
   
   vec = squeeze(v1.img(:,:,slice_no,:));
   p1X = xgrid - vec(:,:,1)/2;
   p1Y = ygrid - vec(:,:,2)/2;
   p2X = xgrid + vec(:,:,1)/2;
   p2Y = ygrid + vec(:,:,2)/2;
   for x=1:v1size(1)
      for y=1:v1size(2)
         line([p1X(x,y) p2X(x,y)], [p1Y(x,y) p2Y(x,y)],'LineWidth', 1.5, 'color', squeeze(abs(vec(x,y,:))));
      end
   end
   
   axis equal
   F = getframe(hfig, [1, 1, figSize, figSize]);
   im = F.cdata;
   imwrite(im, [png_file '.' num2str(slice_no) '.png']);
end

close all 

[xgrid ygrid] = ndgrid(res(1):res(1):v1size(1)*res(1), res(3):res(3):v1size(3)*res(3));
hfig = figure('units','pixels', 'position', [20, 20, figSize, figSize]); %'color', [0 0 0]);
haxes = axes('units','normalized', 'position', [0, 0, 1, 1], 'color', [0 0 0]);
hold on

for slice_no = round(v1size(2)/2-5):round(v1size(2)/2+5)
   bg_img = squeeze(bg.img(:,slice_no,:))';
   imagesc([res(1) v1size(1)*res(1)], [res(3) v1size(3)*res(3)], bg_img, [0 1]);
   colormap(gray)
   set(gca,'ydir','normal');
   hold on;
   
   vec = squeeze(v1.img(:,slice_no,:,:));
   p1X = xgrid - vec(:,:,1)/2;
   p1Y = ygrid - vec(:,:,3)/2;
   p2X = xgrid + vec(:,:,1)/2;
   p2Y = ygrid + vec(:,:,3)/2;
   for x=1:v1size(1)
      for y=1:v1size(3)
         line([p1X(x,y) p2X(x,y)], [p1Y(x,y) p2Y(x,y)],'LineWidth', 1.5, 'color', squeeze(abs(vec(x,y,:))));
      end
   end
   
   axis equal
   F = getframe(hfig, [1, 1, figSize, figSize]);
   im = F.cdata;
   imwrite(im, [png_file '.z.' num2str(slice_no) '.png']);
end




end
