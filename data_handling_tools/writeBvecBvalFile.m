% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function writeBvecBvalFile(bvec, bvec_fname, bval, bval_fname)
% bvec is nDirx3, bval is nDirx1
%
% Usage: 
%    writeBvecBvalFile(bvec, bvec_fname, bval, bval_fname)
%    writeBvecBvalFile(bvec, bvec_fname)

if size(bvec, 2)~=3
   error('bvec should be nDirx3')
end

fid = fopen(bvec_fname, 'w');

for k = 1:size(bvec, 1)
   fprintf(fid, '%22.15f ', bvec(k,:));
   fprintf(fid, '\n');
end
fclose(fid);


if nargin>2
   fid = fopen(bval_fname, 'w');
   fprintf(fid, '%f\n', bval);
   fclose(fid);
end
end
