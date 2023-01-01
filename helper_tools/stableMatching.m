% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [mMatch, wMatch] = stableMatching(mPref, wPref, init_map)
% Men propose version. Implementation of Gale-Shapley "stable marriage" algorithm.
%
% mPref: nxn men's preference. ith row is preference for ith man. Higher value mean higher
%          preference. (i,j) represent the preference-measure of ith man for jth woman.
% wPref: nxn women's preference. jth COLUMN is prerence for jth woman.
% init_map: (optional) nx1 vector of initialization of the  mapping, zero means not mapped. MUST
%           follow the algorithm logic (proposals till the point there is no conflict).
%
% mMatch: nx1 vector of final matches for men. ith man is matched to woman corresponding to ith element.
% wMatch: nx1 vector of final matches for women.
%

% initialy all men and women are single and happy.
mSingle = true(size(mPref,1), 1);
wSingle = true(size(mPref,1), 1);
mMatch = zeros(size(mPref,1), 1);
wMatch = zeros(size(mPref,1), 1);

% Men's proposition matrix - ith row is for ith man. false means not proposed yet, true means already proposed.
propMat = false(size(mPref));

% Fast initialization 
if nargin==3
   mSingle(init_map>0) = false;
   mMatch = init_map;
   wSingle(init_map(init_map>0)) = false;   
   wMatch(init_map(init_map>0)) = find(init_map>0);
   ind = sub2ind(size(propMat), find(init_map>0), init_map(init_map>0));
   propMat(ind) = true;
end

% run while as long as there are free men that hasn't proposed to every woman.
while any(mSingle) && ~all(all(propMat))
   choosenM = find(mSingle, 1); % choose 1st single man, does not matter who proposes first
   
   % choosenM's current best preference, whome he did not propose yet
   bestMatch = max(mPref(choosenM, ~propMat(choosenM,:)));
   MsFavoriteW = find((mPref(choosenM,:)==bestMatch) & ~propMat(choosenM,:), 1);
   propMat(choosenM, MsFavoriteW) = true;
   
   if wSingle(MsFavoriteW) % MsFavoriteW is single, engage them
      mMatch(choosenM) = MsFavoriteW;
      wMatch(MsFavoriteW) = choosenM;
      wSingle(MsFavoriteW) = false;
      mSingle(choosenM) = false;
      
   else % MsFavoriteW is already engaged to someone else
      otherM = wMatch(MsFavoriteW);
      if wPref(otherM, MsFavoriteW) < wPref(choosenM, MsFavoriteW)
         mSingle(choosenM) = false;
         mMatch(choosenM) = MsFavoriteW;
         wMatch(MsFavoriteW) = choosenM;
         
         mSingle(otherM) = true; % otherM is single now
         mMatch(otherM) = 0;
      end
   end
end

% All men and woman are paired
end

