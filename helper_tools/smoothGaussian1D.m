% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function x_smooth = smoothGaussian1D(x, sgma)
% Replacement for smooth() and normcdf() to get rid of some toolbox dependencies. 
%
% Smooths a vector x using convolution with Gaussian pdf. End points are repeated for padding. 
%

xt = 0:6*sgma;
xt = [-1*xt(end:-1:2) xt]; % always odd length

if length(xt)<2
   x_smooth = x;
   
else
   P = normcdfBDP(xt, 0, sgma);
   u = [P(1) diff(P)];
   ulen = length(u)-1; % always even
   temp_x = [ones(ulen,1)*x(1); x(:); ones(ulen,1)*x(end)];
   x_smooth = conv(temp_x, u, 'valid');
   x_smooth = x_smooth(1+ulen/2:length(x)+ulen/2);
end

end

function p = normcdfBDP(x,mu,sigma)
% computes CDF of normal distribution using erfc()
z = (x-mu) ./ sigma;
p = NaN(size(z),class(z));

% Set edge case sigma=0
p(sigma==0 & x<mu) = 0;
p(sigma==0 & x>=mu) = 1;

% Normal cases
if isscalar(sigma)
   if sigma>0
      todo = true(size(z));
   else
      return;
   end
else
   todo = sigma>0;
end
z = z(todo);

% Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))),
% to produce accurate near-zero results for large negative x.
p(todo) = 0.5 * erfc(-z ./ sqrt(2));
end
