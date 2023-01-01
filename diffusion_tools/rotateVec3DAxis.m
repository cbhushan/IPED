% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function vrot = rotateVec3DAxis(v, ax, theta)
% Rotate 3D vectors using Rodrigues' rotation formula
% v - A matrix of 3D vectors, with last dimension of length 3
% ax - A matrix of 3D rotation axis corresponding to each vector in v. ax should have same size as
%      v. Axis in ax need not be of unit length. 
% theta - Angles in radians corresponding to each vector in v. 

if ~isequal(size(v), size(ax))
   error('Size of v and ax must match.')
end
sz = size(v);
if sz(end)~=3
   error('Last dimension of v must have size of 3.')
end

if numel(theta)~=prod(sz(1:end-1))
   error('size of v and theta does not match');
end

ax = reshape(ax, [], 3);
v = reshape(v, [], 3);
theta = theta(:);

% normalize ax
ax_mag = sqrt(sum(ax.^2, 2));
ax = ax./ax_mag(:, [1 1 1]);

costheta = cos(theta);
sintheta = sin(theta);

dotp = sum(ax.*v, 2);
crossp = cross(ax, v, 2);

vrot = v.*costheta(:,[1 1 1]) + crossp.*sintheta(:,[1 1 1]) + (ax.*dotp(:,[1 1 1]).*(1-costheta(:,[1 1 1])));

vrot = reshape(vrot, sz);
clearvars -except vrot

end
