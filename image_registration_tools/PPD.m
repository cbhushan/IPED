% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [W1 W2 W3] = PPD(V1, V2, V3, J1, J2, J3, mask)
%
%  Check PPD_FAST. Its a faster version of this function
%
%
%PPD Rotates the eigenvectors of the diffusion matrix with Preservation of
%Principle component(PPD) algorithm & generates new eigenvectors. Saves the
%output in .nii.gz file.
%
%   V1 - Principal eigenvector (corresponding to eigenvalue L1)
%   V2 - Second eigenvector (corresponding to eigenvalue L2)
%   V3 - Third eigenvector (corresponding to eigenvalue L3)
%       such that L1>L2>L3
%
%   J1, J2, J3 - Jacobian of each component (x,y,z) of map
%   mask - 3D mask for operation (optional)
%
%   W1 - Rotated principal eigenvector
%   W2 - Rotated second eigenvector
%   W3 - Rotated third eigenvector
%

warning('Calling PPD_FAST. Its faster version of this function.');
[W1 W2 W3] = PPD_FAST(V1, V2, V3, J1, J2, J3);
return;
   
fprintf('\n==== Applying PPD ====\n');

% convert to double for higher precision
V1.img = double(V1.img);
V2.img = double(V2.img);
V3.img = double(V3.img);

J1 = double(J1);
J2 = double(J2);
J3 = double(J3);

[len bre dep dim4] = size(V1.img);

if exist('mask')==0
	mask = ones(len, bre, dep);
end

W1 = struct('hdr',V1.hdr, 'filetype',V1.filetype, 'machine',V1.machine, 'original',V1.original);
W1.fileprefix = [V1.fileprefix '.PPD'];
W1.img = zeros(len, bre, dep, dim4);

W2 = struct('hdr',V2.hdr, 'filetype',V2.filetype, 'machine',V2.machine, 'original',V2.original);
W2.fileprefix = [V2.fileprefix  '.PPD'];
W2.img = zeros(len, bre, dep, dim4);

W3 = struct('hdr',V3.hdr, 'filetype',V3.filetype, 'machine',V3.machine, 'original',V3.original);
W3.fileprefix = [V3.fileprefix '.PPD'];
W3.img = zeros(len, bre, dep, dim4);

for m = 1:len
	for n = 1:bre
		for o = 1:dep
			
			switch mask(m,n,o)
				case 1
					f1 = squeeze(J1(m, n, o, :)); % 1st row of Jacobian
					f2 = squeeze(J2(m, n, o, :)); % 2nd row of Jacobian
					f3 = squeeze(J3(m, n, o, :)); % 3rd row of Jacobian
					F = [f1' ; f2' ; f3'] + eye([3,3]);
					
					e1 = squeeze(V1.img(m, n, o, :));
					e2 = squeeze(V2.img(m, n, o, :));
					e3 = squeeze(V3.img(m, n, o, :));
					
					e1 = e1 / norm(e1);
					e2 = e2 / norm(e2);
					e3 = e3 / norm(e3);
					
					tempVec = F*e1;
					tempVecNorm = norm(tempVec);
					switch tempVecNorm
						case 0     % check for zero vector
							n1 = zeros(size(tempVec));
						otherwise
							n1 = (tempVec)/tempVecNorm;  % convert to unit vector
					end
					
					tempDotProduct = dot(e1, n1);
					
					% Eliminate the posibility of complex theta2
					switch (tempDotProduct > 0)
						case 1
							cosTheta1 = min(tempDotProduct, 1.0);
						case 0
							cosTheta1 = max(tempDotProduct, -1.0);
					end
					
					axis1 = cross(e1,n1);  % e1 to n1; not unit vector
					R1 = vrrotvec2mat_cos([axis1' cosTheta1]); 
					
					tempVec = F*e2;
					tempVecNorm = norm(tempVec);
					switch tempVecNorm
						case 0     % check for zero vector
							n2 = zeros(size(tempVec));
						otherwise
							n2 = (tempVec)/tempVecNorm;  % convert to unit vector
					end
					
					% Projection of n2 perpendicular to n1
					Proj_n2 = n2 - dot(n2,n1)*n1;
					Proj_n2 = (Proj_n2)/norm(Proj_n2);
					
					e2RotR1 = R1*e2;
					axis2 = cross(e2RotR1, Proj_n2); % R1*e2 to Proj_n2; not unit vector
					tempDotProduct = dot(Proj_n2, e2RotR1);
					
					% Eliminate the posibility of complex theta2
					switch (tempDotProduct > 0)
						case 1
							cosTheta2 = min(tempDotProduct, 1.0);
						case 0
							cosTheta2 = max(tempDotProduct, -1.0);
					end
					
					R2 = vrrotvec2mat_cos([axis2' cosTheta2]);
					
					R = R2*R1;

					W1.img(m,n,o,:) = R*e1;
					W2.img(m,n,o,:) = R*e2;
					W3.img(m,n,o,:) = R*e3;
			end
		end
	end
	fprintf('PPD: %2.2f %%\n',m/len*100)
end

save_nii_gz(W1, []);
save_nii_gz(W2, []);
save_nii_gz(W3, []);

fprintf('Done PPD..\n');
end

