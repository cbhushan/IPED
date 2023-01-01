% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function lbl_sorted = sortLabelsNearestCentroids(lbl, x, ref_lbl, ref_x)
% Sorts (actually maps) input labels to make them nearest (euclidean sense) to the labels in
% reference.
%
% lbl - nx1 labels
% x - nxd coordinates
% ref_lbl - mx1 reference labels
% ref_x - mxd coordinates

sz = size(x);

lbl_u = unique(lbl, 'stable');
n_lbl = length(lbl_u);

ref_lbl_u = unique(ref_lbl, 'stable');
n_ref_lbl = length(ref_lbl_u);


% find centroids
lbl_centroid = zeros(n_lbl, sz(2));
for k = 1:n_lbl
   lbl_centroid(k, :) = mean(x(lbl==lbl_u(k), :), 1);
end

ref_lbl_centroid = zeros(n_ref_lbl, sz(2));
for k = 1:n_ref_lbl
   ref_lbl_centroid(k, :) = mean(ref_x(ref_lbl==ref_lbl_u(k), :), 1);
end

D = pdist2(lbl_centroid, ref_lbl_centroid, 'euclidean');
[~, target_ind] = min(D, [], 2);

% make label mapping unique, if required
[C, ia] = unique(target_ind, 'stable');
if length(C)~=length(target_ind)
   if n_lbl>n_ref_lbl
      error('Not implemented yet! May be swap ref and target.')
   else
      unassigned_target_inds = setdiff(1:n_ref_lbl, C);
      repeated_ind = setdiff(1:n_lbl, ia);
      target_ind(repeated_ind) = unassigned_target_inds(1:length(repeated_ind));      
   end
end

lbl_sorted = zeros(size(lbl));
for k = 1:n_lbl
   lbl_sorted(lbl==lbl_u(k)) = ref_lbl_u(target_ind(k));
end

end
