% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function tf = isRAS_sform(s, res)
% Checks if sform matrix is in cannonical form or not.
%   s - A nifti struct or 4x4 matrix
%   res - 3x1 vector representing resolution, ignored when s is struct

if isstruct(s)
   if s.hdr.hist.sform_code == 0
      error('sform_code is unset! Try using add_sform first.')
   end
   
   sformT = zeros([3 4]);
   sformT(1,:) = s.hdr.hist.srow_x;
   sformT(2,:) = s.hdr.hist.srow_y;
   sformT(3,:) = s.hdr.hist.srow_z;
   res = abs(s.hdr.dime.pixdim(2:4));
   M = sformT(1:3, 1:3) * diag(1./res);
   
else
   M = s(1:3, 1:3) * diag(1./res(1:3));
end


reorient_matrix = zeros(3,3);
B = sort(abs(M(:)), 'descend');
max_mask = abs(M)>=B(3);
temp = M(max_mask);
reorient_matrix(max_mask) = temp./abs(temp);

tf = isequal(reorient_matrix, eye(3));

end
