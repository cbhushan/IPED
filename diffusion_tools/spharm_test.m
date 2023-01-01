% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function spharm_test( )
% using spherical harmonics table

clc

N = 100;
phiRadians = linspace(-pi/2, pi/2, N);
thetaRadians = linspace(0, 2*pi, N);
[thetaRadians, phiRadians] = meshgrid(thetaRadians, phiRadians);

shDegree = 0;
shOrder = 0;
Y = spharm(shDegree, shOrder, thetaRadians, phiRadians);
Ycomp = repmat(1/(2*sqrt(pi)), [N N]);
err = Y-Ycomp;
max(abs(err(:)))


clear Ycomp Y err 
shDegree = 1;
shOrder = -1;
Y = spharm(shDegree, shOrder, thetaRadians, phiRadians);
Ycomp = sqrt(3/(8*pi))*(exp(-1i*thetaRadians).*cos(phiRadians));
err = Y-Ycomp;
max(abs(err(:)))


clear Ycomp Y err
shDegree = 1;
shOrder = 0;
Y = spharm(shDegree, shOrder, thetaRadians, phiRadians);
Ycomp = sqrt(3/(4*pi))*sin(phiRadians);
err = Y-Ycomp;
max(abs(err(:)))


clear Ycomp Y err
shDegree = 1;
shOrder = 1;
Y = spharm(shDegree, shOrder, thetaRadians, phiRadians);
Ycomp = -1*sqrt(3/(8*pi))*(exp(1i*thetaRadians).*cos(phiRadians));
err = Y-Ycomp;
max(abs(err(:)))

clear Ycomp Y err
shDegree = 2;
shOrder = -2;
Y = spharm(shDegree, shOrder, thetaRadians, phiRadians);
Ycomp = sqrt(15/(32*pi))*(exp(-2i*thetaRadians).*(cos(phiRadians)).^2);
err = Y-Ycomp;
max(abs(err(:)))

clear Ycomp Y err
shDegree = 2;
shOrder = 1;
Y = spharm(shDegree, shOrder, thetaRadians, phiRadians);
Ycomp = -1*sqrt(15/(8*pi))*(exp(1i*thetaRadians).*(cos(phiRadians)).*sin(phiRadians));
err = Y-Ycomp;
max(abs(err(:)))


clear Ycomp Y err
shDegree = 3;
shOrder = -2;
Y = spharm(shDegree, shOrder, thetaRadians, phiRadians);
Ycomp = sqrt(105/(32*pi))*(exp(-2i*thetaRadians).*(cos(phiRadians)).^2.*sin(phiRadians));
err = Y-Ycomp;
max(abs(err(:)))


clear Ycomp Y err
shDegree = 6;
shOrder = -2;
Y = spharm(shDegree, shOrder, thetaRadians, phiRadians);
Ycomp = sqrt(1365/pi)/64*(exp(-2i*thetaRadians).*(cos(phiRadians)).^2.*(33*sin(phiRadians).^4-18*sin(phiRadians).^2+1));
err = Y-Ycomp;
max(abs(err(:)))

end

