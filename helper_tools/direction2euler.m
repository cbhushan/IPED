% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [thetaD, rot_mat] = direction2euler(pos)
% pos is Nx1 vector - typically be positions returned by generateIPEDdirSingleShell()

epsilon = 1e-12; % Minimum value to treat a number as zero
nDirections = size(pos, 1);

% xvec = [1 0 0]
xvec = zeros(nDirections, 3);
xvec(:,1) = 1;

axis_rot = cross(xvec, pos, 2);
rot_mat = zeros(3, 3, nDirections);

% remove numerical issues
pos(abs(pos)<epsilon) = 0;
axis_rot(abs(axis_rot)<epsilon) = 0;

% find rotation matrices
theta = zeros(nDirections, 3);
for k = 1:nDirections
   r = [axis_rot(k,:) pos(k,1)];
   rot_mat(:,:,k) = vrrotvec2mat_cos(r);
   theta(k,:) = rotmat2euler(rot_mat(:,:,k));
end

% in degrees
thetaD = theta*(180/pi);

% take care of [-1 0 0]
ind = find(abs(pos(:,1)+1)<1e-5);
thetaD(ind,:) = 0;
thetaD(ind,3) = 180;

end


function theta = rotmat2euler(R)
% convertes 3x3 rot_mat to three euler angles
% based on R = Rx * Ry * Rz - as in par2affineMat.m

theta = zeros(1,3);
theta(1) = atan2(-1*R(2,3), R(3,3));
theta(2) = atan2(R(1,3), sqrt(R(2,3)^2 + R(3,3)^2));
theta(3) = atan2(-1*R(1,2), R(1,1));

end


% function theta = decompose_rotation(R)
% % based on R = Rz * Ry * Rx;
% 
% x = atan2(R(3,2), R(3,3));
% y = atan2(-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)));
% z = atan2(R(2,1), R(1,1));
% theta = [x,y,z];
% end

