% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [M, unM, Minv, unMinv] = createMaskOperators(msk)
% Returns masking and unmasking operators for input mask. The operators should be applied on indexed
% vectors.
%   M - Masks data 
%   unM - Unmasks the masked data. Fills in zeros outside the mask. 
%   Minv - Masks data using inverse of the input mask. 
%   unMinv - Unmasks the data masked by Minv. Fills in zeros outside the mask. 
%

msk = msk~=0;
sz = size(msk);

M = speye(prod(sz));
M(~msk, :) = [];

if nargout>1
   unM = speye(prod(sz));
   unM(:, ~msk) = [];
end

if nargout>2
   Minv = speye(prod(sz));
   Minv(msk, :) = [];
end

if nargout>3
   unMinv = speye(prod(sz));
   unMinv(:, msk) = [];
end

end
