% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [Iout, shft_best, MI] = max_mutual_information_circshift(Is, Is_mask, Im, Im_mask, circshift_sizes)
% compute MI for all possible combinations & return the one with max MI

if ~isequal(size(Is), size(Im))
   size(Is)
   size(Im)
   error('Size of these images do not match!')
end

if ~isequal(ndims(Is), size(circshift_sizes, 2))
   error('Number of dimensions does not match length of circshift_vectors.')
end

% histeq to make it easier to find max
Im_histeq = normalize_intensity(Im, [0 100], Im_mask>0);
Is = normalize_intensity(Is, [0 100], Is_mask>0);
Is(Is_mask>0) = histeqSmooth(Is(Is_mask>0), Im_histeq(Im_mask>0), 256, 24);

MI = zeros(size(circshift_sizes, 1),1);
for ishift = 1:size(circshift_sizes, 1)
   shft = circshift_sizes(ishift, :);
   
   img_temp = circshift(Im_histeq, shft);
   msk_temp = circshift(Im_mask, shft);
   
   msk = msk_temp>0 & Is_mask>0;
   MI(ishift) = mutualInfoParzen(img_temp(msk), Is(msk), 250, [0 1], 12, 8);
end

% find max
ind = find(MI==max(MI));
shft_best = circshift_sizes(ind(1), :);

% make return image
Iout = circshift(Im, shft_best);

end

