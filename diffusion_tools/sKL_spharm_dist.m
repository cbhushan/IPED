% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function sKL_dist = sKL_spharm_dist(S1, bMatrices1, S2, bMatrices2, shOrder)
% computes sKL distance (Chiang 2007) between two diffusion-weighted images
% S1 & S2. Last dimension is the no of diffusion directions. Expects to
% include b=0 image as the first image in S1 & S2
%

S1Size = size(S1);
S2Size = size(S2);

if ~isequal(S1Size(1:end-1), S2Size(1:end-1))
   error('D1 & D2 should have same no of voxels')
end
nVoxels = prod(S2Size(1:end-1));

S1 = double(S1);
S2 = double(S2);

% throw away negative numbers
S1(S1<0) =  0;
S2(S2<0) =  0;


% Find spharm coeff for D1
c1 = fit_ADC(S1, bMatrices1, shOrder);

temp = reshape(bMatrices1, 9, []);
temp = sum(temp.^2, 1);
b0mask = (temp==0);
clear temp

S1 = reshape(S1, nVoxels, []);
c1 = reshape(c1, nVoxels, []);

% extract b=0 data
if sum(b0mask)==1
   S1_0 = S1(:, b0mask);
elseif sum(b0mask)>1
   S1_0 = mean(S1(:, b0mask), 2);
end
S1(:, b0mask) = [];
bMatrices1(:,:,b0mask) = [];
clear b0mask

% Find log of ADC
nDirection = size(bMatrices1, 3);
D1log = zeros(size(S1));
for iDir = 1:nDirection
   [~, bval] = sorted_eig(bMatrices1(:,:,iDir));
   D1log(:, iDir) = log( abs(log(S1(:,iDir)./S1_0))/bval(1,1) );
end
d1 = fit_data(D1log, bMatrices1, shOrder);
clear S1 S1_0 bMatrices1 D1log nDirection


% Find spharm coeff for D2
c2 = fit_ADC(S2, bMatrices2, shOrder);

temp = reshape(bMatrices2, 9, []);
temp = sum(temp.^2, 1);
b0mask = (temp==0);
clear temp

S2 = reshape(S2, nVoxels, []);
c2 = reshape(c2, nVoxels, []);

% extract b=0 data
if sum(b0mask)==1
   S2_0 = S2(:, b0mask);
elseif sum(b0mask)>1
   S2_0 = mean(S2(:, b0mask), 2);
end
S2(:, b0mask) = [];
bMatrices2(:,:,b0mask) = [];
clear b0mask

% Find log of ADC
nDirection = size(bMatrices2, 3);
D2log = zeros(size(S2));
for iDir = 1:nDirection
   [~, bval] = sorted_eig(bMatrices2(:,:,iDir));
   D2log(:, iDir) = log( abs(log(S2(:,iDir)./S2_0))/bval(1,1) );
end
d2 = fit_data(D2log, bMatrices2, shOrder);
clear S2 S2_0 bMatrices2 D2log nDirection


% compute sKL distance
term1 = sum(c1.*(d1-d2), 2) ./ c1(:,1);
term2 = sum(c2.*(d2-d1), 2) ./ c2(:,1);

sKL_dist = (term1+term2)/2;


if numel(S1Size(1:end-1))==1
   S1Size(end) = 1; % add 1 to make it reshape-able
   S1Size(end+1) = 1;
end
sKL_dist = reshape(sKL_dist, S1Size(1:end-1));

end

