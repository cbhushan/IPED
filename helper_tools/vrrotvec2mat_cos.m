% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [ m ] = vrrotvec2mat_cos ( r )
%VRROTVEC2MAT_COS Similar to vrrotvec2mat_dti(described below), with
%following difference:
%
%   Last dimension of row vector r is NOT angle. It is cosine of the angle.
%   This might be useful to avoid numerical errors associated with acos
%   following by cos.
%
%   VRROTVEC2MAT_DTI Convert rotation from r-angle to matrix representation
%
%   m = vrrotvec2mat(r) returns a matrix representation of the rotation
%   defined by the axis-cosine(angle) rotationvector, r.
%
%   The rotation vector, r, is a row vector of four elements, where the
%   first three elements specify the rotation axis (may NOT be unit
%   vector), and the last element defines the angle (in radians).
%
%   Note: Simulink 3D Animation toolbox has another function (vrrotvec2mat)
%   with same functionality, but possibly diffrent implementation. This function uses
%   matrix form of Rodrigues' rotation formula.

axis = r(1:3);

if all(axis == 0)
	warning('Zero vector as rotation axis. Will generate identity matrix as rotation matrix.');
	m = eye(3);
else
	axis = axis/norm(axis);
	tensorProduct = axis' * axis;
	crossProdMat = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
	
	%rotation matrix
	cosAngle = r(4);
	sinAngle = sqrt(1-(cosAngle^2));
	m = tensorProduct + cosAngle*(eye(size(tensorProduct))-tensorProduct) + sinAngle*crossProdMat;
end

end

