% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [epi_scl, bfield, bfield_msk, alpha] = matchLocalIntensity(epiImg, sImg, epi_mask, struct_mask, opts)
% This function tries to remove bias field effect in b=0 image by comparing local histograms to
% MPRAGE image. Number of local sub-divisions is defined by opts.ndiv. 
%
% When opts.method=1 then each local subdivision in epiImg is adjusted to have closest to identity
%    intensity mapping to sImg.  
% When opts.method=2 then each local subdivision in epiImg is adjusted to have closest to the global
%    mapping to sImg. 
%
% bfield is the smooth version of scaling computed for each sub-division. 
% alpha is the original scaling computed for each sub-division and each intensity level
%


defaultoptions = struct(...
   'nbins', 200, ...
   'parzen_width', 20, ...
   'ndiv', [2 3 3], ...
   'method', 1);

if(~exist('opts','var')),
   opts = defaultoptions;
else
   tags = fieldnames(defaultoptions);
   for i=1:length(tags)
      if(~isfield(opts,tags{i})),  opts.(tags{i})=defaultoptions.(tags{i}); end
   end
   if(length(tags)~=length(fieldnames(opts))),
      error('Unknown options found.');
   end
end

if ~isequal(size(epiImg), size(sImg))
   error('size of epiImg and sImg must be same');
end

nbins = opts.nbins;
parzen_width = opts.parzen_width;
ndivX = opts.ndiv(1);
ndivY = opts.ndiv(2);
ndivZ = opts.ndiv(3);

% stip volume tightly
sz_in = size(sImg);
[tb1, tb2, tb3] = find_bounding_box(epi_mask>0|struct_mask>0, 0);
epiImg = epiImg(tb1,tb2,tb3);
sImg = sImg(tb1,tb2,tb3);
epimsk = epi_mask(tb1,tb2,tb3)>0;
strmsk = struct_mask(tb1,tb2,tb3)>0;
tsz = size(sImg);
tmsk = false(tsz);

% Refine local histogram matching
xind = floor(linspace(1,tsz(1)+1, ndivX+1));
yind = floor(linspace(1,tsz(2)+1, ndivY+1));
zind = floor(linspace(1,tsz(3)+1, ndivZ+1));

[~,~,c0,b0] = histcParzen(1-epiImg(epimsk), [0 1], nbins, parzen_width); 
m = b0>=0 & b0<=1;
c0 = c0(m)/sum(c0);

% compute inverse mapping
[~, hMap0inv_pp] = histeqSmooth(sImg(strmsk), 1-epiImg(epimsk), nbins, parzen_width);

alpha = zeros(ndivX,ndivY,ndivZ, length(c0));
bscl = zeros(tsz);
for xx = 1:ndivX
   for yy = 1:ndivY
      for zz = 1:ndivZ
         m1 = tmsk;
         m1(xind(xx):xind(xx+1)-1, yind(yy):yind(yy+1)-1, zind(zz):zind(zz+1)-1) = true;
         
         if any(epimsk(:) & m1(:) & strmsk(:))
            [~,~,c1,b1] = histcParzen(1-epiImg(epimsk & m1), [0 1], nbins, parzen_width);
            c1 = c1/sum(c1);
            
            [~, ~, hMap1, t1] = histeqSmooth(1-epiImg(epimsk & m1), sImg(strmsk & m1), nbins, parzen_width);
            t1 = t1(m);
            hMap1 = hMap1(m);
            c1 = c1(m);
            
            t1 = t1(end:-1:1);
            if opts.method == 1
               alpha(xx,yy,zz,:) = ((1-hMap1)./t1).*c1;
            else
               alpha(xx,yy,zz,:) = ((1-ppval(hMap0inv_pp, hMap1))./t1).*c1;
            end
         else
            alpha(xx,yy,zz,:) = c0;
         end
         temp = alpha(xx,yy,zz,:);
         bscl(m1) = sum(temp(isfinite(temp)));
      end
   end
end

% smooth the scale map
alpha(~isfinite(alpha)) = 0;
bsclS = imgaussian(bscl, 4);
bsclS = imgaussian(bsclS, 4);
bsclS = imgaussian(bsclS, 4);

epi_scl = epiImg.*bsclS;
epi_scl = epi_scl/max(epi_scl(epimsk));
bfield = epi_scl./epiImg;
bfield_msk = isfinite(bfield);
bfield(~bfield_msk) = 1;

% change output size to original size
epi_scl = resize_to_orig(epi_scl, sz_in, tb1, tb2, tb3);
bfield = resize_to_orig(bfield, sz_in, tb1, tb2, tb3);
bfield_msk = resize_to_orig(bfield_msk, sz_in, tb1, tb2, tb3);

end

function vol_out = resize_to_orig(vol, out_sz, tb1, tb2, tb3)
vol_out = zeros(out_sz);
vol_out(tb1, tb2, tb3) = vol;
end


