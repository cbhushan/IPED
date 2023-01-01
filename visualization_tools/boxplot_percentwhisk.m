% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function bh = boxplot_percentwhisk(data, percen, varargin)
% For compact boxplot use: bh = boxplot_percentwhisk(data, percen, Tg, G, ..., 'plotstyle', 'compact')
%
% Allows data to be a cell vector with each cell having different number of elements. It is
% converted to a 2D matrix by repeating each cell as required. See cellDataRepeat() below. 
% NOTE: This uses LCM to repeat the data - can result into huge data matrix!
%
% for each grouping along the x-axis, replot the whisker to be based on
% percentiles of the data. Also, only the outliers outside the new whiskers
% are plotted. For the percentile, it is calculated and then rounded up to
% the next nearest data value. For example, if the 98th percentile is
% between the 4th highest and 5th highest data value, then the percentile
% is set to be the 4th highest data value.

%keyboard

data = cellDataRepeat(data);

if nargin>2
   tf = cellfun(@isstr, varargin);
   compact_boxplot = ismember('compact', varargin(tf));
else
   compact_boxplot = false;
end

if compact_boxplot
   
   % first plot out the Matlab version of boxplot so that you can get the
   % handles.
   bh = boxplot(data, varargin{2:end});
   Tg = varargin{1};
   
   % get the necessary handles.
   b = get(gca,'children');
   c = get(b,'children');
   UW = findobj(c,'Tag','Whisker');
   sz = length(Tg);
   
   for k = 1:sz
      msk = ismember(varargin{2}, Tg(sz-k+1)); % somehow index in UW is reverse of the original order      
      v = percentile(data(msk), percen);
      
      data_lowpercent(k) = v(1);
      data_highpercent(k) = v(2);
      UW_lims = get(UW(k),'ydata');
      set(UW(k),'ydata',[data_lowpercent(k) data_highpercent(k)])      
   end
   
else % normal boxplot
   
   % first plot out the Matlab version of boxplot so that you can get the
   % handles.
   bh = boxplot(data, varargin{:});
   
   
   % get the necessary handles.
   b = get(gca,'children');
   c = get(b,'children');
   UW = findobj(c,'Tag','Upper Whisker');
   UAJ = findobj(c,'Tag','Upper Adjacent Value');
   OUT = findobj(c,'Tag','Outliers');
   LW = findobj(c,'Tag','Lower Whisker');
   LAJ = findobj(c,'Tag','Lower Adjacent Value');
   
   for i=1:size(data,2)
      v = percentile(data(:,i), percen);
      data_lowpercent(i) = v(1);
      data_highpercent(i) = v(2);
      UW_lims = get(UW(size(data,2)-i+1),'ydata');
      set(UW(size(data,2)-i+1),'ydata',[UW_lims(1) data_highpercent(i)])
      set(UAJ(size(data,2)-i+1),'ydata',[data_highpercent(i) data_highpercent(i)])
      LW_lims = get(LW(size(data,2)-i+1),'ydata');
      set(LW(size(data,2)-i+1),'ydata',[data_lowpercent(i) LW_lims(2)])
      set(LAJ(size(data,2)-i+1),'ydata',[data_lowpercent(i) data_lowpercent(i)])
      %     outliers = get(OUT(size(data,2)-i+1),'ydata');
      %     outlier_index = find(outliers == data_highpercent(i));
      %     outliers(1:outlier_index) = data_highpercent(i);
      %     set(OUT(size(data,2)-i+1),'ydata',outliers);
   end
end

rng = max(data_highpercent)-min(data_lowpercent);
temp = ylim();
ylim([0 max(data_highpercent)+0.1*rng]);

end


function [v] = percentile(data, kpercent)

y = sort(data, 'ascend');
k = kpercent./100;
ind = round(k * length(y));

ind(ind<1) = 1;
ind(ind>length(data)) = length(data);
v = y(ind);

end


function d = cellDataRepeat(cdata)
% Uses LCM to repeat the data - can result into huge data matrix!

if iscell(cdata)   
    ncols = numel(cdata);
   if ncols==1
      d = cdata{1};
      return;
   end
   
   sz = cellfun(@numel, cdata);
   L = lcm(sz(1), sz(2));
   for k = 3:ncols
      L = lcm(L, sz(k));
   end
   
   d = zeros(L, ncols);
   for k = 1:ncols
      n_rep = L/sz(k);
      
      temp = cdata{k};
      temp = temp(:);
      temp = temp(:, ones(1, n_rep));      
      d(:,k) = temp(:);
   end   
   
else % regular input
   d = cdata;
end

end

