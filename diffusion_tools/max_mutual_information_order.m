% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [Iout, ord_best, MI] = max_mutual_information_order(Is, Is_mask, Im, Im_mask, img_order)
% compute MI for all possible combinations & return the one with max MI

if ~isequal(sort(size(Is)), sort(size(Im)))
   size(Is)
   size(Im)
   error('Size of these images can never match irrespective of any ordering!')
end

% histeq to make it easier to find max MI
Im_histeq = normalize_intensity(Im, [0 100], Im_mask>0);
Is = normalize_intensity(Is, [0 100], Is_mask>0);
Is(Is_mask>0) = histeqSmooth(Is(Is_mask>0), Im_histeq(Im_mask>0), 256, 24);

MI = zeros(size(img_order, 1),1);
for iorder = 1:size(img_order, 1)
   ord = img_order(iorder, :);
   
   img_temp = Im_histeq;
   msk_temp = Im_mask;
   
   if ord(1)<0
      img_temp = flipdim(img_temp, 1);
      msk_temp = flipdim(msk_temp, 1);
   end
   
   if ord(2)<0
      img_temp = flipdim(img_temp, 2);
      msk_temp = flipdim(msk_temp, 2);
   end   
   
   if ord(3)<0
      img_temp = flipdim(img_temp, 3);
      msk_temp = flipdim(msk_temp, 3);
   end
   
   img_temp = permute(img_temp, abs(ord));
   msk_temp = permute(msk_temp, abs(ord));
   
   
   if isequal(size(img_temp), size(Is))
      msk = msk_temp>0 & Is_mask>0;
      MI(iorder) = mutualInfoParzen(img_temp(msk), Is(msk), 250, [0 1], 16, 8);
   else
      MI(iorder) = -0.1;
   end
   
end

% find max
ind = find(MI==max(MI));
ord_best = img_order(ind(1), :);

% make return image
img_temp = Im;

if ord_best(1)<0
   img_temp = flipdim(img_temp, 1);
end

if ord_best(2)<0
   img_temp = flipdim(img_temp, 2);
end

if ord_best(3)<0
   img_temp = flipdim(img_temp, 3);
end

Iout = permute(img_temp, abs(ord_best));


end

