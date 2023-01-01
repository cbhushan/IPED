% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function GFA = GFA_dwi(dwi)
% Computes GFA from DWIs. First fit DWIs with SH basis and then uses coefficients to compute GFA.

dwiSize = size(dwi);
n_vox = prod(dwiSize(1:end-1));

dwi = reshape(dwi, n_vox, []);

GFA = std(dwi, 0, 2) ./ rms(dwi,2);
GFA(~isfinite(GFA)) = 0;

GFA = reshape(GFA, dwiSize(1:end-1));

end
