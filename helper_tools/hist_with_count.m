% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [ bins, counts ] = hist_with_count( data, nBins)
% Plots histogram with count on top of each bin.

if ~exist('nBins','var')
  nBins = 10;
end

% get the edges, bin centers
edges = linspace(min(data),max(data),nBins+1); %# edges go from minimum to maximum of distribution
bins = (edges(1:end-1)+edges(2:end))/2;

% get the counts and the bin-index
[counts,binIdx] = histc(data,edges);
counts(end-1) = sum(counts(end-1:end));  %# combine last two bins
counts(end) = [];                        %# 
binIdx(binIdx==nBins+1) = nBins;         %# also fix the last bin index

% plot the counts and bins (not edges) with `bar`
bar(bins,counts);

% Set the axes limits such that you have enough space for the labels
ylim([0,2*max(counts)]);

% add the labels. Vertically align such that the text goes from the y-coordinate
% down (as opposed to being centered on the y-coordinate).
for b = 1:nBins
    text(bins(b)-19,counts(b)+98, num2str(counts(b)) ,'VerticalAlignment','top', 'Units', 'data')
end

end

