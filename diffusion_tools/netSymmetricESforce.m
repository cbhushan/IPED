% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [F, F_tangent] = netSymmetricESforce(pos)
% Computes net electrostatic repulsion force at each position on sphere due to all other points and
% phantom points located at diametrically opposite location of each point. 
%
% POS is nx3 matrix where each row represents a point on unit sphere/n-ball/circle - make sure that
% each row has unit L2 norm!
%

[n, m] = size(pos);
[a,b] = ndgrid(1:n, 1:n); % b(:) represents nth row of pos & is in sorted order
C = [b(:) a(:)]; % All combination + repeation
C(b(:)==a(:),:) = []; % throw away repeatations

d1 = pos(C(:,1),:) - pos(C(:,2),:); % distance 1
d1_m = sum(d1.^2, 2).^(3/2);
f1 = d1./d1_m(:, ones(1, m));

d2 = pos(C(:,1),:) + pos(C(:,2),:); % distance from diametrically opposite point
d2_m = sum(d2.^2, 2).^(3/2);
f2 = d2./d2_m(:, ones(1, m));

% Sum all forces on each pos; See the hack below
a = (n-1):(n-1):(n*(n-1));
temp = cumsum(f1, 1);
f1_sum = diff([zeros(1, m); temp(a,:)], 1, 1); % because b(:) & C(:,1) is in sorted order! 

temp = cumsum(f2, 1);
f2_sum = diff([zeros(1, m); temp(a,:)], 1, 1);

F = f1_sum + f2_sum;

if nargout>1
   dotp = sum(F.*pos, 2);
   F_tangent = F - pos.*dotp(:, ones(1, m));
end

end

% A beautiful summation hack
%
% If A = [3,4,3] and B is vector of length(sum(A)) then following will produce a vector C where
% C(1) is the sum of the first A(1) elements of B, C(2) is the sum of the next A(2) elements of B, etc. 
% foo = cumsum(B);
% C = diff([0 foo(cumsum(A))])
