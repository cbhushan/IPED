% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [lbl_sorted, lbl_sorted_non_unique, pref_msr_all, map_unique, map_non_unique] = matchLabelsMaxPrefMeasure(lbl, ref_lbl, pref_msr)
% Maps input labels to ref-labels such that it maximizes the preference measure defined by input
% pref_msr. Assumes that lbl and ref_lbl have same size. When number of unique labels in ref_lbl is
% less than lbl, then lbl may be mapped to some dummy labels which do not exist in ref_lbl.
%
% Uses Gale-Shapley "stable marriage" algorithm, to resolve conflicts between same choice (lbl
% prpopses). This also means following:
%   * lbl gets largest possible preference measure among all possible solutions of the mapping.
%   * From ref_lbl point of view, it gets worst possible preference measure among all possible
%     solutions of the mapping.
%   * Swapping lbl and ref_lbl can result into different mapping.
%
% Inputs:
%    lbl - nx1 vector of labels
%    ref_lbl - nx1 vector reference labels
%    pref_msr - string defining label-wise preference measures. Supported measures are 'dice' and 'cardinality'.
%
% Outputs:
%    lbl_sorted - nx1 vector of mapped labels. Each label will be have one-to-one mapping.
%    lbl_sorted_non_unique - nx1 vector of mapped labels. It can have many-to-one mappings.
%    pref_msr_all - Preference matrix used for exhaustive search.
%    map_unique - mx2 matrix of mapping used for lbl_sorted. map_unique(:,1) is unique labels in lbl. 
%                 map_unique(:,2) is the map corresponding to map_unique(:,1).
%    map_non_unique - same as map_unique but for lbl_sorted_non_unique.
%

lbl = lbl(:);
ref_lbl = ref_lbl(:);
sz = size(lbl);
pref_msr = lower(pref_msr);

if ~isequal(sz, size(ref_lbl))
   error('lbl and ref_lbl must have same number of elements.')
end

if ~ismember(pref_msr, {'dice', 'cardinality'})
   error('pref_msr must be either of ''dice'' or ''cardinality''.')
end

lbl_u = unique(lbl, 'stable');
n_lbl = length(lbl_u);

ref_lbl_u = unique(ref_lbl, 'stable');
n_ref_lbl = length(ref_lbl_u);

if n_lbl>n_ref_lbl
   n_dummy_lbls = n_lbl - n_ref_lbl;
   dummy_ref_label = [1:n_dummy_lbls] + max(ref_lbl_u);
   ref_lbl_u(end+1:end+n_dummy_lbls) = dummy_ref_label;
   n_ref_lbl = length(ref_lbl_u);   
   warning('matchLabelsMaxPrefMeasure:NotEnoughRefLabels', ...
      '%d dummy labels added for the mapping.', n_dummy_lbls);
end

% find pref measure
pref_msr_all = zeros(n_lbl, n_ref_lbl);
for i = 1:n_lbl
   msk_i = (lbl==lbl_u(i));
   Ni = sum(msk_i);
   for r = 1:n_ref_lbl      
      msk_r = (ref_lbl==ref_lbl_u(r));
      switch pref_msr
         case 'dice'
            pref_msr_all(i, r) = 2*sum(msk_i & msk_r)/(Ni + sum(msk_r));
         case 'cardinality'
            pref_msr_all(i, r) = sum(msk_i & msk_r);
      end
   end
end

[~, target_ind] = max(pref_msr_all, [], 2);
target_ind_non_unique = target_ind;

% make label mapping unique, if required
[C, ia] = unique(target_ind, 'stable');
if length(C)~=length(target_ind)
   init_map = zeros(size(target_ind));
   init_map(ia) = C;
   target_ind = stableMatchingSolution(pref_msr_all, init_map);
end

lbl_sorted = zeros(size(lbl));
for i = 1:n_lbl
   lbl_sorted(lbl==lbl_u(i)) = ref_lbl_u(target_ind(i));
end

lbl_sorted_non_unique = zeros(size(lbl));
for i = 1:n_lbl
   lbl_sorted_non_unique(lbl==lbl_u(i)) = ref_lbl_u(target_ind_non_unique(i));
end

% return the mapping
map_unique = [lbl_u(:) ref_lbl_u(target_ind(:))];
map_non_unique = [lbl_u(:) ref_lbl_u(target_ind_non_unique(:))];

end


function target_ind = stableMatchingSolution(pref_msr_all, init_map)
% convert pref_msr_all to a square matrix and then use stable matching

sz = size(pref_msr_all);

if sz(1)>sz(2)
   error('This should have been sorted out earlier! Something went wrong!')
   
elseif sz(1)<sz(2)
   more_rows = zeros(sz(2)-sz(1), sz(2));
   pref_msr_all = cat(1, pref_msr_all, more_rows);
   init_map = cat(1, init_map(:), zeros(sz(2)-sz(1), 1));
end

target_ind = stableMatching(pref_msr_all, pref_msr_all, init_map);
target_ind = target_ind(1:sz(1));
end


