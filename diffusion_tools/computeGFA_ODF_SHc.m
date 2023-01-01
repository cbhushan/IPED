% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function GFA = computeGFA_ODF_SHc(c, Pf)
% Computes GFA from (modified) SH coefficients of ODFs.
% c  - first dim of c should be SH coefficients. SHc should be ordered 'naturally'. Other dim are spatial.
%
% Pf - scaling parameter for (SH coefficients of DWI-fit); Eg: 2*pi*P_k(0) term in eq(1) of the paper
%    below. Vector of same length as dim-1 of c. 
%
% GFA - computed GFA values. Length of first dim would be always 1. 
%
% J. Cohen-Adad et al., "Detection of multiple pathways in the spinal cord using q-ball imaging."
%      NeuroImage, 42 (2008): 739 - 749, DOI: http://dx.doi.org/10.1016/j.neuroimage.2008.04.243  
%

cSize = size(c);
c = reshape(c, cSize(1), []);
[nterms, nvox] = size(c);

% Normalize for scaling parameter - get SHc for DWI data fit
Pf = Pf(:);
c = c ./ Pf(:, ones(1, nvox));

% GFA
c = c.^2;
GFA = sqrt(1 - (c(1,:)./sum(c,1)));

GFA(~isfinite(GFA)) = 0;
GFA = reshape(GFA, [1 cSize(2:end)]);

end
