% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [D,Dp] = createDWithPeriodicBoundary(N1,N2)
% [D,Dp] = createD(N1,N2)
%
% Generates the sparse two-dimensional finite difference (first-order neighborhoods) matrix D 
% for an image of dimensions N1xN2 (rows x columns).  The optional output
% argument Dp is the transpose of D.  Also works for vector images if one
% of N1 or N2 is 1.
%
% Justin Haldar (jhaldar@usc.edu) 11/02/2006
% Inspired by the dgrad function, written by W. C. Karl 4/98

if (not(isreal(N1)&&(N1>0)&&not(N1-floor(N1))&&isreal(N2)&&(N2>0)&&not(N2-floor(N2))))
    error('Inputs must be real positive integers');
end
if ((N1==1)&&(N2==1))
    error('Finite difference matrix can''t be generated for a single-pixel image');
end

D1 = [];
D2 = [];

if (N1 > 1)&&(N2>1)    
    e = ones(N1,1);
    if (numel(e)>2)
        T = spdiags([e,-e],[0,1],N1,N1);
        T(N1,1)=-1;
        E = speye(N2);
        D1 = kron(E,T);
    end
    e = ones(N2,1);
    if (numel(e)>2)
        T = spdiags([e,-e],[0,1],N2,N2);
        T(N2,1)=-1;
        E = speye(N1);
        D2 = kron(T,E);
    end
else % Image is a vector of length max(N1,N2)
    N = max(N1,N2);
    e = ones(N,1);
    D1 = spdiags([e,-e],[0,1],N,N);
    D1(end,1) = -1;
end

D = [D1;D2];

if (nargout > 1)
    Dp = D';
end

end
