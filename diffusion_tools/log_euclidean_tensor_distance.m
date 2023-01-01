% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [ d ] = log_euclidean_tensor_distance(eig_file1, eig_file2)
% Returns Log-Euclidean (tensor) distance between DTI images. eig_file1 & eig_file2 are
% *.eig.nii.gz files as saved by estimate_tensor.m  
%
% If A & B are two tensors then Log-Euclidean distance is defined as 
% (Arsigny et al. MRM 2006) 
%
% dist(A,B) = sqrt(trace((logm(A)-logm(B))^2))
%
% Eigen decomposition of diffusion tensor (Symmetric positive definite) is exploited to
% simplify computation.
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

% compute log 
e1(:,10:12) = log(abs(e1(:,10:12)));
e2(:,10:12) = log(abs(e2(:,10:12)));

% logm(A)-logm(B)
e_diff = eig2tensor(e1) - eig2tensor(e2);

% compute sqrt( trace(e_diff^2) )
term =  sqrt( ...
          e_diff(:,1).^2 + 2*(e_diff(:,2).^2) + 2*(e_diff(:,3).^2) + ...
          e_diff(:,4).^2 + 2*(e_diff(:,5).^2) +    e_diff(:,6).^2 ...
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

