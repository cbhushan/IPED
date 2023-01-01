% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [ diagN ] = n_shore(radial_order )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

F = radial_order / 2;
n_c = round(1 / 6.0 * (F + 1) * (F + 2) * (4 * F + 3));
diagN = zeros(n_c,1);

counter = 1;
for n = 0:radial_order,
    for l = 0:2:n,
        for m = -l:l,
            diagN(counter) = (n * (n + 1))^2;
            counter = counter + 1;
        end;
    end;
end;
diagN=diag(diagN);

end

