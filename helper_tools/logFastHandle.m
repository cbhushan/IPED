% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function h = logFastHandle(bins)
% Returns a handle to compute fast logarithm based on lookup tables.
% bins - number of bins to be used for joint-histogram

tbl = load('log_lookup.mat');

% for range1
l = length(tbl.log_lookup1);
md1_rng1 = ones(bins*bins,1)*((l-1)/diff(tbl.range1));
md2_rng1 = ones(2*bins, 1)*((l-1)/diff(tbl.range1));
md1_const1 = ones(bins*bins,1) - (tbl.range1(1)*(l-1)/diff(tbl.range1));
md2_const1 = ones(2*bins,1) - (tbl.range1(1)*(l-1)/diff(tbl.range1));

% for range2
l = length(tbl.log_lookup2);
md1_rng2 = ones(bins*bins,1)*((l-1)/diff(tbl.range2));
md2_rng2 = ones(2*bins, 1)*((l-1)/diff(tbl.range2));
md1_const2 = ones(bins*bins,1) - (tbl.range2(1)*(l-1)/diff(tbl.range2));
md2_const2 = ones(2*bins,1) - (tbl.range2(1)*(l-1)/diff(tbl.range2));

h = @logFast;

   function log_data = logFast(data, mode)
      % mode = 1, when data is bins x bins (for joint histogram)
      % mode = 2, when data is 2*bins (for concatenated marginal histogram)
      
      log_data = zeros(size(data));
      
      ind_0 = (data<tbl.range1(1));
      ind_1 = (data<tbl.range1(2)) & ~ind_0;
      ind_2 = (data<=tbl.range2(2)) & ~(ind_1 | ind_0);
      ind_3 = (data>tbl.range2(2));
      
      log_data(ind_0) = log(tbl.range1(1)/2); % constant for lowest quadrant
      
      if mode == 1
         ind = uint32((data(ind_1).*md1_rng1(ind_1)) + md1_const1(ind_1));
         log_data(ind_1) = tbl.log_lookup1(ind);
         
         ind = uint32((data(ind_2).*md1_rng2(ind_2)) + md1_const2(ind_2));
         log_data(ind_2) = tbl.log_lookup2(ind);
         
      else % mode = 2
         ind = uint32((data(ind_1).*md2_rng1(ind_1)) + md2_const1(ind_1));
         log_data(ind_1) = tbl.log_lookup1(ind);
         
         ind = uint32((data(ind_2).*md2_rng2(ind_2)) + md2_const2(ind_2));
         log_data(ind_2) = tbl.log_lookup2(ind);
      end
      
      log_data(ind_3) = log(data(ind_3)); % outside lookup table
   end
end

