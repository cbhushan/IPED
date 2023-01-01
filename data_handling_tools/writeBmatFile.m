% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function writeBmatFile(bMatrices, fname)
% bMatrices should be 3x3xnDir

fid = fopen(fname, 'w');
for iDir = 1:size(bMatrices,3)
   temp = squeeze(bMatrices(:,:,iDir))';
   temp = temp(:)';
   fprintf(fid, '%22.15f %22.15f %22.15f\n%22.15f %22.15f %22.15f\n%22.15f %22.15f %22.15f\n\n', temp);
end
fclose(fid);
end
