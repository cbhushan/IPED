% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [mi, nmi] = mutualInfoParzen(x, y, nbins, range, window_size, nthreads, adjust_LH)
% Computes mutual information between two vectors x, y. x & y must be of equal size. range defines
% the [low high] range of the intensity values in x & y. First computes joint histogram (parzen
% estimate) of nbins x nbins. window_size is parzen window size in units of bins, which can be
% either one or two element vector (when two elements, 1st element defines the window_size for x
% and 2nd element for y). nthreads is number of parallel threads to be used for histogram 
% computation.
%
%  Usage: 
%       [mi, nmi] = mutualInfoParzen(x, y, nbins)
%       [mi, nmi] = mutualInfoParzen(x, y, nbins, range) % range = [low high]
%       [mi, nmi] = mutualInfoParzen(x, y, nbins, range, window_size)
%       [mi, nmi] = mutualInfoParzen(x, y, nbins, range, window_size, nthreads)
%       [mi, nmi] = mutualInfoParzen(x, y, nbins, range, window_size, nthreads, adjust_LH)
%
% Also see test_mutualInfoParzen()

if nargin == 3
   window_size = max(2, round(double(nbins)/10));
   nthreads = 12;
   adjust_LH = true;
   range = [0 1];
   
elseif nargin == 4
   window_size = max(2, round(double(nbins)/10));
   nthreads = 12;
   adjust_LH = true; 
   
elseif nargin == 5
   nthreads = 12;
   adjust_LH = true;
   
elseif nargin == 6
   adjust_LH = true;
end

if numel(window_size)==1
   window_size = window_size * [1 1];
end

if adjust_LH
   l = 0 - max(window_size)/nbins/2;
   h = 1 + max(window_size)/nbins/2;
else
   l = 0;
   h = 1;
end

x = clipData(x, range(1), range(2));
y = clipData(y, range(1), range(2));

hist12 = mutual_histogram_parzen_two_variable_size_multithread_double(double(x(:)), double(y(:)), l, h, ...
   double(nbins), double(window_size(1)), double(window_size(2)), double(nthreads));

hist1 = sum(hist12, 1);
hist2 = sum(hist12, 2);

% log computation can be slow. May be use lookup table.
c = log(eps);
log_hist12 = ones(size(hist12))*c;
log_hist1 = ones(size(hist1))*c;
log_hist2 = ones(size(hist2))*c;

ind = find(hist1 > eps);
log_hist1(ind) = log(hist1(ind));

ind = find(hist2 > eps);
log_hist2(ind) = log(hist2(ind));

ind = find(hist12 > eps);
log_hist12(ind) = log(hist12(ind));


h1 = hist1(:) .* log_hist1(:);
h2 = hist2(:) .* log_hist2(:);
h12 = hist12 .* log_hist12;

mi = sum(h12(:)) - sum(h1) - sum(h2);

if nargout>1
   nmi = (sum(h1) + sum(h2)) / sum(h12(:)); % studholme et al. 1999
end

end


function I1out = clipData(I1, low, high)
   I1out = double((I1-low)/(high-low));
   I1out(I1out<0) = 0;
   I1out(I1out>1) = 1;
end
