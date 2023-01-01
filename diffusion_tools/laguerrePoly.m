% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


% Calculates Generalized Laguerre polynomial
% Szegö: Orthogonal Polynomials, 1958, (5.1.10)
function L = laguerrePoly(n,a,x)
l=zeros(n+1);          
if(n==0)
    l(1,:)=1;
else
    l(1,:)=[zeros(1,n), 1];
    l(2,:)=[zeros(1, n-1), -1, (a+1)];
    for i=3:n+1
        A1 = 1/(i-1) * (conv([zeros(1, n-1), -1, (2*(i-1)+a-1)], l(i-1,:)));
        A2 = 1/(i-1) * (conv([zeros(1, n), ((i-1)+a-1)], l(i-2,:)));
        B1=A1(length(A1)-n:1:length(A1));
        B2=A2(length(A2)-n:1:length(A2));
        l(i,:)=B1-B2;
     end;
end
L=polyval(l(n+1,:),x); 
