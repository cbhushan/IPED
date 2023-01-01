% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [hist12, h, cbar_axes] = scatterDensity(d1, d2, d1range, d2range, nbins, winsz, nCPU, logstr)
% Usage: 
%     scatterDensity(d1, d2)
%     scatterDensity(d1, d2, d1range, d2range, nbins, winsz, nCPU, logstr)
%     [hist12, h, cbar_axes] = scatterDensity(...)


rng = [-0.1 1.1];

if numel(d1) ~= numel(d2)
   error('Number of elements in d1 and d2 must match.')
end

if ~exist('nbins', 'var')
   nbins = 256;
end

if ~exist('winsz', 'var')
   winsz = 12;
end

if ~exist('nCPU', 'var')
   nCPU = 10;
end

if exist('d1range', 'var')
   d1 = double((d1(:)-d1range(1))/(d1range(2)-d1range(1)));
else
   [d1, l, h] = normalize_intensity(d1(:), [0 100]);
   d1range = [l h];
end

if exist('d2range', 'var')
   d2 = double((d2(:)-d2range(1))/(d2range(2)-d2range(1)));
else
   [d2, l, h] = normalize_intensity(d2(:), [0 100]);
   d2range = [l h];
end

hist12 = mutual_histogram_parzen_variable_size_multithread_double(d1(:), d2(:), rng(1), rng(2), nbins, winsz, nCPU);
hist12 = hist12*numel(d1);

rngV = linspace(rng(1), rng(2), nbins);
b1 = (rngV*(d1range(2)-d1range(1))) + d1range(1);
b2 = (rngV*(d2range(2)-d2range(1))) + d2range(1);

figure;
if exist('logstr', 'var') && strcmpi(logstr, 'log')
   h = imagesc(b1, b2,  log(hist12));
else
   h = imagesc(b1, b2,  hist12);
end
set(gca,'YDir','normal');
xlabel('d1');
ylabel('d2'); 
% grid on; 

cbar_axes = colorbar;
v = caxis;
if exist('logstr', 'var') && strcmpi(logstr, 'log')
   v(1) = log(0.2); % counts should be usually more than 0.1
   caxis(v);
   vsp = linspace(v(1), v(2), 5);
   vspstr = cell(length(vsp),1);
   for cc = 1:length(vsp)
      vspstr{cc} = sprintf('%d',round(exp(vsp(cc))));
   end
   
else
   vsp = linspace(v(1), v(2), 5);
   vspstr = cell(length(vsp),1);
   for cc = 1:length(vsp)
      vspstr{cc} = sprintf('%d',vsp(cc));
   end
end

set(cbar_axes, 'YTick', vsp, 'YTickLabel', vspstr);

if nargout<1
   clear hist12
end

end
