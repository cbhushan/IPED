% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [ d ] = sKL_tensor_distance(eig_file1, eig_file2)
% Returns symmetric KL distance between DTI images. Same as sqrt of J-divegence. 
% eig_file1 & eig_file2 are *.eig.nii.gz files as saved by estimate_tensor.m  
%
%     T1, T2: Must be symmetric positive definite matrix (3x3)
%             These are covariance matrix of gaussian random
%             variables (say p & q respectively).
%
%     Then, J = (KL(p|q) + KL(q|p))/2, and
%             = 0.5*sqrt( trace(T1inv*T2 + T2inv*T1) - 2*n )
%           distance = sqrt(J)
%           
%     This function uses the closed form obtained after simiplifying the
%     expressions. (Yoshizawa and Tanabe, SUT J. Math 1999; Wang and Vemuri, TMI 2005)
%

if ischar(eig_file1)
   eig1 = load_untouch_eig_gz(eig_file1);
elseif ~isfield(eig_file1,'untouch') || eig_file1.untouch ~= 1
   error('Please use ''load_untouch_eig_gz.m'' to load the structure.');
else
   eig1 = eig_file1;
   clear eig_file1;
end

if ischar(eig_file2)
   eig2 = load_untouch_eig_gz(eig_file2);
elseif ~isfield(eig_file2,'untouch') || eig_file2.untouch ~= 1
   error('Please use ''load_untouch_eig_gz.m'' to load the structure.');
else
   eig2 = eig_file2;
   clear eig_file2;
end



% sanity check of sizes 
if ~isequal(size(eig1.img), size(eig2.img))
   error('Dimension size does not match in eig_file1 & eig_file2');
end
sz = size(eig1.img);
mask = (eig1.img(:,:,:,10)>0 & eig1.img(:,:,:,11)>0 & eig1.img(:,:,:,12)>0) | ...
       (eig2.img(:,:,:,10)>0 & eig2.img(:,:,:,11)>0 & eig2.img(:,:,:,12)>0);

mask4D = mask(:,:,:, ones(1,12));

e1 = double(reshape(eig1.img(mask4D>0), sum(mask(:)), 12));
e2 = double(reshape(eig2.img(mask4D>0), sum(mask(:)), 12));

% make sure eigenvalues are +ve
e1(:,10:12) = abs(e1(:,10:12));
e2(:,10:12) = abs(e2(:,10:12));

% compute inverse
e1inv = e1;
e2inv = e2;
e1inv(:,10:12) = 1./e1(:,10:12);
e2inv(:,10:12) = 1./e2(:,10:12);

T1 = eig2tensor(e1); clear e1;
T2 = eig2tensor(e2); clear e2;
T1inv = eig2tensor(e1inv); clear e1inv;
T2inv = eig2tensor(e2inv); clear e1inv;

% compute sqrt( trace(T1inv*T2 + T2inv*T1)-6 )
term =  0.5.*sqrt( ...
          T1inv(:,1).*T2(:,1) + 2*(T1inv(:,2).*T2(:,2)) + 2*(T1inv(:,3).*T2(:,3)) + ...
          T1inv(:,4).*T2(:,4) + 2*(T1inv(:,5).*T2(:,5)) +    T1inv(:,6).*T2(:,6)  + ...
          ...
          T2inv(:,1).*T1(:,1) + 2*(T2inv(:,2).*T1(:,2)) + 2*(T2inv(:,3).*T1(:,3)) + ...
          T2inv(:,4).*T1(:,4) + 2*(T2inv(:,5).*T1(:,5)) +    T2inv(:,6).*T1(:,6) ...
          - 6 ...
         );

d = zeros(sz(1:3));
d(mask>0) = term;

end




function T = eig2tensor(eig)
% eigen decomposition to symmetric notation of matrix
%  
%  | 1 2 3 |
%  |   4 5 |
%  |     6 |
%

T(:,1) = eig(:,10).*(eig(:,1).^2) + ...
         eig(:,11).*(eig(:,4).^2) + ...
         eig(:,12).*(eig(:,7).^2);

T(:,2) = eig(:,10) .* eig(:,1) .* eig(:,2) + ...
         eig(:,11) .* eig(:,4) .* eig(:,5) + ...
         eig(:,12) .* eig(:,7) .* eig(:,8);

T(:,3) = eig(:,10) .* eig(:,1) .* eig(:,3) + ...
         eig(:,11) .* eig(:,4) .* eig(:,6) + ...
         eig(:,12) .* eig(:,7) .* eig(:,9);

T(:,4) = eig(:,10).*(eig(:,2).^2) +...
         eig(:,11).*(eig(:,5).^2) + ...
         eig(:,12).*(eig(:,8).^2);

T(:,5) = eig(:,10) .* eig(:,2) .* eig(:,3) + ...
         eig(:,11) .* eig(:,5) .* eig(:,6) + ...
         eig(:,12) .* eig(:,8) .* eig(:,9);

T(:,6) = eig(:,10).*(eig(:,3).^2) +...
         eig(:,11).*(eig(:,6).^2) + ...
         eig(:,12).*(eig(:,9).^2);
end

