% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function filt_vol = anisotropic_diffusion_filter3d(imgvol, niter, options)
% Performs anisotropic diffusion filtering of input 3D image volume. For consistent results
% intensities in imgvol should be normalized and resolution information should be used. 
%
% Arguments:
%   imgvol - input image volume
%   niter - number of iteration
%   options - optional input 
%        options.kappa  - diffusion coefficient
%        options.lambda - size of time step; 0<lambda<1/3
%        option.eq - 1; Perona, Malik diffusion equation No 1 (favours high contrast edges over low
%                       contrast ones)
%                    2; Perona, Malik diffusion equation No 2 (favours wide regions over smaller
%                       ones)
%
% Reference:
%  * P. Perona and J. Malik. Scale-space and edge detection using ansotropic diffusion
% IEEE Transactions on Pattern Analysis and Machine Intelligence,
% 12(7):629-639, July 1990.
%
%  * Mackiewich, B. Intracranial Boundary Detection and Radio Frequency Correction in Magnetic
% Resonance Images; Simon Fraser University, 1995 (masters thesis)
%
%  * DW.Shattuck et al. Magnetic Resonance Image Tissue Classification Using a Partial Volume Model 
%    NeuroImage, 2001, 13, 856-876
%
%

opt = struct( ...
   'kappa', 25, ...
   'lambda', 1/7,...
   'eq', 1, ...
   'voxRes', [1 1 1] ...
   );

if ndims(imgvol)~=3
   error('imgvol must be a 3D volume');
end

% set options using options
if exist('options', 'var')
   fnames = fieldnames(opt);
   for iname = 1:length(fnames)
      if isfield(options, fnames{iname})
         opt = setfield(opt, fnames{iname}, getfield(options, fnames{iname}));
      end
   end
end

deltaX = 1/opt.voxRes(1);
deltaY = 1/opt.voxRes(2);
deltaZ = 1/opt.voxRes(3);

imgvol = double(imgvol);
filt_vol = imgvol;
kappa = opt.kappa;

for iters = 1:niter
   
   dE = convn(filt_vol, [0 -deltaX deltaX], 'full'); 
   dE = dE(:,2:size(dE,2)-1,:);
   
   dW = convn(filt_vol, [deltaX -deltaX 0], 'full'); 
   dW = dW(:,2:size(dW,2)-1,:);
   
   dN = convn(filt_vol, [0; -deltaY; deltaY], 'full'); 
   dN = dN(2:size(dN,1)-1,:,:);
   
   dS = convn(filt_vol, [deltaY; -deltaY; 0], 'full'); 
   dS = dS(2:size(dS,1)-1,:,:);
   
   kernel = zeros(1,1,3); kernel(2) = -deltaZ; kernel(3) = deltaZ;
   dU = convn(filt_vol, kernel, 'full'); 
   dU = dU(:,:,2:size(dU,3)-1);
   
   kernel = zeros(1,1,3); kernel(1) = deltaZ; kernel(2) = -deltaZ;
   dD = convn(filt_vol, kernel, 'full'); 
   dD = dD(:,:,2:size(dD,3)-1);
   
   if opt.eq == 1
      filt_vol = filt_vol + opt.lambda * (...
         (exp(-1*(dE/kappa).^2) .* dE) + (exp(-1*(dW/kappa).^2) .* dW) + ...
         (exp(-1*(dN/kappa).^2) .* dN) + (exp(-1*(dS/kappa).^2) .* dS) + ...
         (exp(-1*(dU/kappa).^2) .* dU) + (exp(-1*(dD/kappa).^2) .* dD));
      
   elseif opt.eq == 2
      filt_vol = filt_vol + opt.lambda * (...
         (dE./(1+(dE/kappa).^2)) + (dW./(1+(dW/kappa).^2)) + ...
         (dN./(1+(dN/kappa).^2)) + (dS./(1+(dS/kappa).^2)) + ...
         (dU./(1+(dU/kappa).^2)) + (dD./(1+(dD/kappa).^2)));
   else
      error('Invalid input for diffusion equation. It must be either 1 or 2.')
   end
end

end
