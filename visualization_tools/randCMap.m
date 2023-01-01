% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function cout = randCMap(sz)
% Randomly generate colormap
h = 0:1/6:5/6;

if sz<=8
   h = [1/3 h(1:sz-2)];
   s = [0 ones(1,sz-2)];
   v = [1   ones(1,sz-2)];
   
elseif sz<=14
   h = [1/3 h h];
   h = h(1:sz-1);
   s = [0 ones(1,sz-2)];
   v = [1   ones(1,6) 0.5*ones(1,6)];
   v = v(1:sz-1);
   
elseif sz<=20
   h = [1/3 h h h];
   h = h(1:sz-1);
   s = [0 ones(1,12) 0.5*ones(1,6)];
   s = s(1:sz-1);
   v = [1   ones(1,6) 0.5*ones(1,6) ones(1,6)];
   v = v(1:sz-1);
   
else
   h = [1/3 h h h];
   s = [0  ones(1,12) 0.5*ones(1,6)];
   v = [1   ones(1,6) 0.5*ones(1,6) ones(1,6)];
   
   h = [h rand(1, sz-20)];
   s = [s (rand(1, sz-20)*0.7)+0.3];
   v = [v (rand(1, sz-20)*0.7)+0.3];
   
end
c = hsv2rgb([h(:) s(:) v(:)]);
p = randperm(size(c,1));
cout = cat(1, [0.1 0.1 0.1], c(p,:));

end

