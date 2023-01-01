% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [V, D] = sorted_eig(A)
% Returns sorted (descending) eigenvalue and eigenvector of matrix A

[V1, D1] = eig(A);

D = diag(sort(diag(D1),'descend'));
[~, ind] = sort(diag(D1),'descend');
V = V1(:,ind);

end
