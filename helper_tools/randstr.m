% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function rand_str = randstr(n)
% generates a random string of lower case letters and numbers of length n


if usejava('jvm')
   tmp_name = strrep(char(java.util.UUID.randomUUID),'-','');
   while length(tmp_name)<n
      tmp_name = [tmp_name strrep(char(java.util.UUID.randomUUID),'-','')];
   end
else
   tmp_name = num2str(feature('timing','cpucount'));
   while length(tmp_name)<n
      tmp_name = [tmp_name num2str(feature('timing','cpucount'))];
   end
end

rand_str = tmp_name(1:n);

end
