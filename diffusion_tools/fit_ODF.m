% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function sh_c = fit_ODF(dwi, bMatrices, shOrder, user_options)
% Estimates coefficients for SH basis functions (real & symmetric) for ODF estimated after 
% Funk-Radon Transform (FRT)
% 
% dwi - Diffusion weighted data. Last dimension is the no of diffusion encoding directions(n).
%       eg: X*Y*Z*n or V*n
%
% bMatrices - 3*3*n
% shOrder - SH order to fit (even positive integer)
%
% sh_c - coefficients of first nTerms of modified SH. Last dimention is the
%        estimated coefficients. Eg: X*Y*Z*nTerms or V*nTerms
%

opt = struct( ...
   'xi', 0.34, ...
   'odf_lambda', 0.006, ...
   'odf_type', 'FRT'... % FRT / FRACT
   );

% set options using user_options
if exist('user_options', 'var')
   fnames = fieldnames(opt);
   for iname = 1:length(fnames)
      if isfield(user_options, fnames{iname})
         opt = setfield(opt, fnames{iname}, getfield(user_options, fnames{iname}));
      end
   end
end


dwi = double(dwi);
bMatrices = double(bMatrices);

dwiSize = size(dwi);
if isvector(dwi)
   s0Size = 1;
   nVoxels = 1;
else
   s0Size = dwiSize(1:end-1);
   nVoxels = prod(s0Size);
end

if dwiSize(end)~=size(bMatrices, 3)
   error('Number of diffusion images in dwi and size of bMatrices does not match.')
end

dwi = reshape(dwi, nVoxels, []);
nDir = dwiSize(end);

% compute bvals
bval = zeros(1,nDir);
for iDir = 1:nDir
   bval(iDir) = trace(bMatrices(:,:,iDir));
end

% Find diffusion encoding direction (same as largest eigen vector of bMatrices)
ind = find(bval);
q = zeros(numel(ind),3);
for i = 1:numel(ind)
    [V,D] = eig(bMatrices(:,:,ind(i)));
    [~,tmp] = max(diag(D));
    q(i,:) = V(:,tmp);
end


% FRT and FRACT Computation
HarmonicOrder = shOrder;
xi = opt.xi;
lambda = opt.odf_lambda;

[S,L] = sph_harm_basis([q(:,1),q(:,2),q(:,3)],HarmonicOrder,2);

for order = 0:2:HarmonicOrder
    lpp = legendre(order,xi);
    lpo = legendre(order,0);
    lpn = legendre(order,-xi);
    Pg(order+1) = 2*lpo(1)-lpp(1)-lpn(1);
    Pf(order+1) = 2*lpo(1);
end
sphericalHarmonicMatrixReg = (S'*S+lambda*diag(L.^2.*(L+1).^2))\(S');

% throw away b=0 images
dwi = dwi(:,ind);
dwi = permute(dwi, [2,1]);
outSize = [size(S,2) nVoxels];

if strcmpi(opt.odf_type, 'FRACT')
   fractMatrixReg = diag(2*pi*Pg(L+1))*sphericalHarmonicMatrixReg;
   sh_c = fractMatrixReg * dwi;
   
elseif strcmpi(opt.odf_type, 'FRT')
   frtMatrixReg = diag(2*pi*Pf(L+1))*sphericalHarmonicMatrixReg;
   sh_c = frtMatrixReg * dwi;
   
else
   error('Incorrect option for odf_type')
end

sh_c = permute(sh_c, [2,1]);
sh_c = reshape(sh_c, s0Size, []);

end


% OLD FRT implementation 
% % modified SH basis
% nDirection = size(dwi, 1);
% vec = zeros(nDirection, 3);
% for iDir = 1:nDirection
%    [V, ~] = sorted_eig(bMatrices(:,:,iDir));
%    vec(iDir,:) = V(:,1);
% end
% [theta, phi, ~] = cart2sph(vec(:,1),vec(:,2),vec(:,3));
% [Y, nTerms, L] = spharm_Hardi(shOrder, theta, phi); % basis matrix
% 
% Pf =  zeros(shOrder+1,1);
% for order = 0:2:shOrder
%    lpo = legendre(order,0);
%    Pf(order+1) = 2*lpo(1);
% end
% 
% LB = diag(L.^2 .* (L+1).^2); % Laplace Beltrami smoothing matrix
% FRT = diag(2*pi*Pf(L+1)); % FRT matrix
% 
% % Old implementation - correct but slow
% % r = 1;
% % for order = 0:2:shOrder
% %    for m = -order:order
% %       LB(r,r) = order^2*(order+1)^2;
% %       FRT(r,r) = 2*pi*(-1)^(order/2)*factorial(order)/(4^(order/2)*factorial(order/2)^2);
% %       r = r+1;
% %    end
% % end
% 
% c = FRT * ((Y'*Y + lambda*LB)\Y') * dwi;
% c = reshape(c', [s0Size nTerms]);



