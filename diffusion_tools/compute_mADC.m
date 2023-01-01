% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function adc = compute_mADC(dwi, DEout)
% dwi - dxn 2D matrix; d is number of diffusion weighting, n is number of voxels
% DEout - Output of checkDiffusionEncodingScheme()
% adc - mean ADC, 1xn

dwi = double(abs(dwi));
b0_mean = mean(dwi(DEout.zero_bval_mask,:), 1);

% throw away b=0 images
dwi = dwi(~DEout.zero_bval_mask,:);
bval = DEout.bval(~DEout.zero_bval_mask);

nDir = size(dwi, 1);
adc = zeros(size(b0_mean));
dwi_c = zeros(size(b0_mean)); % non_zero counter
for k = 1:nDir
   temp = abs(-1*log(dwi(k,:)./b0_mean)/bval(k));   
   msk = isfinite(temp);   
      
   % use only dwis which are >0 for stable estimate (& keep track of it)
   temp(~msk) = 0;
   dwi_c = dwi_c + double(msk);   
   adc = adc + temp;
end
adc(~isfinite(adc)) = 0;
adc = abs(adc)./dwi_c; % mean ADC

end
