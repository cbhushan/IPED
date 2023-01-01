% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function h = errorbar_percentile(Xin, data, line_type, prct, varargin)
% wrapper around errorbar() to enable customized location of error bar and line plots. Based on the
% inputs (Xin, data, line_type, prct) it is converted to errorbar(X,Y,L,U, ...). 1st dim (or
% columns) of data represents the several overvation for a Xin. That is, 2nd dim of data should be
% same as first-dim (or length in case of vector) of Xin. 3rd dim of data represents data for
% another line in errorbar plots. 
% 
% Xin - vector or 2D matrix for x-axis of plot
% data - 2D or 3D matrix of data. 
%        When 2D : Xin must be vector. Second dim of data should be same as length of Xin.
%        When 3D : Xin can be 1D or 2D. Second dim of data should be same as first dim of Xin. When
%                  Xin is 2D, 3rd dim of data must be of same size as 2nd dim of Xin. Xin is
%                  replicated approrpriately when it is 1D.
%
% line_type - String, must of either of 'mean', 'median' or 'mode'. Always computed along 1st dim of
%             data. 
% prct - [low high] percentile for error bars. Always computed along 1st dim of data. 
%
% varargin - other errorbar optional input
%

if nargin<4
   error('Must have atleast 4 inputs.')
end

if isvector(Xin)
   X = Xin(:);
end
szX = size(X);


switch line_type
   case 'mean'
      Y = squeeze(mean(data, 1));
   case 'median'
      Y = squeeze(median(data, 1));
   case 'mode'
      Y = squeeze(mode(data, 1));
   otherwise
      error('Unrecognized line_type: %s', line_type)
end
szY = size(Y);

if szX(1)~=szY(1)
   error('1st dimension of X and 2nd dim of data are not appropriate. See usage.')
elseif szX(2)==1
   X = X(:, ones(1,szY(2)));
elseif szX(2)~=szY(2)
   error('2nd dimension of X and 3rd dim of data are not appropriate. See usage.')
end


if numel(prct)~=2
   error('prct must be vector of length 2')
end
E = prctile(data, prct, 1);
L = squeeze(E(1,:,:));
U = squeeze(E(2,:,:));

h = errorbar_percentile(X, Y, L, U, varargin);

end


