% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [D,Dp] = createAnisoDwithMask(N1,N2,N3,msk,anisoWt,resWt)
% Generates sparse matrix operator D for an image of dimensions N1xN2xN3 such that it computes
% weighted finite difference only using voxel inside the input mask i.e. boundary voxels (of mask or
% volume) for which the "forward" voxel lies outside mask/volume is not used.
%
% Similar to createDwithMask, but allows to specify different weights at each voxel.
% Applied weight for D is the mean of weights at corresponding voxels.
%
% resWt is a weight vector of length 3, which represents weight due to voxel dimensions. Optional
% input - when not present same weight is applied to all dimension.
%
% Dp is the transpose of D.

if (not(isreal(N1)&&(N1>0)&&not(N1-floor(N1))&&isreal(N2)&&(N2>0)&&not(N2-floor(N2))))
    error('Inputs must be real positive integers');
end

if ((N1==1)&&(N2==1)&&(N3==1))
    error('Finite difference matrix can''t be generated for a single-pixel image');
end

if ~isequal(size(msk), [N1, N2, N3])
   error('Mask size does not match!')
end

if ~isequal(size(anisoWt), [N1, N2, N3])
   error('size of anisotropic wt. does not match!')
end

if ~exist('resWt', 'var')
   resWt = [1 1 1];
end

msk = msk>0; 
msk_ind = find(msk); clear msk;
m = N1*N2*N3;

% D along dim 1
ind_incr = 1;
edge_mask = false(N1,N2,N3);
edge_mask(end,:,:) = true;
D = DinMask(msk_ind, edge_mask, ind_incr, m, anisoWt);

% D along dim 2
ind_incr = N1;
edge_mask = false(N1,N2,N3);
edge_mask(:,end,:) = true;
D2 = DinMask(msk_ind, edge_mask, ind_incr, m, anisoWt);
D = [D.*resWt(1); D2.*resWt(2)];
clear D2;

% D along dim 3
ind_incr = N1*N2;
edge_mask = false(N1,N2,N3);
edge_mask(:,:,end) = true;
D3 = DinMask(msk_ind, edge_mask, ind_incr, m, anisoWt);
clear anisoWt edge_mask msk_ind

D = [D; D3.*resWt(3)];
clear D3;

% remove tonnes of zeros rows
temp = any(D,2);
D = D(temp,:);

if (nargout > 1)
    Dp = D';
end

clearvars -except D Dp
end

function D = DinMask(msk_ind, edge_mask, ind_incr, m, wt)
edge_ind = find(edge_mask);
i = setdiff(msk_ind, edge_ind); % throw away volume edge voxels
clear edge_ind

Lia = ismember(i+ind_incr, msk_ind);
i = i(Lia);

s = (wt(i) + wt(i+ind_incr))/2; % mean of weights at corresponding voxels
j = [i; i+ind_incr];
i = [i; i];
s = [s; -s];
D = sparse(i,j,s,m,m);
clearvars -except D
end
