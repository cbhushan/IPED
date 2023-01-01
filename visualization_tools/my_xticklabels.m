% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function ht = my_xticklabels(varargin)
%
% modified from my_yticklabels to work with xticklabels
%
%MY_YTICKLABELS replaces YTickLabels with "normal" texts
%   accepting multiline texts and TEX interpreting
%   and shrinks the axis to fit the texts in the window
%
%    ht = my_xticklabels(Ha, xtickpos, xtickstring)
% or
%    ht = my_xticklabels(xtickpos, xtickstring)
%
%  in:    xtickpos     YTick positions [N*1]
%        xtickstring   Strings to use as labels {N*1} cell of cells
%
% Examples:
% plot(randn(20,1))
% ytl = {{'one';'two';'three'} '\alpha' {'\beta';'\gamma'}};
% h = my_yticklabels(gca,[-1 0 1],ytl);
% % vertical
% h = my_yticklabels([-1 0 1],ytl, ...
%     'Rotation',-90, ...
%     'VerticalAlignment','middle', ...
%     'HorizontalAlignment','center');

% Pekka Kumpulainen 5.3.2008

textopts = {};
if length(varargin{1})==1 && ...
        ishandle(varargin{1}) && ...
        strcmpi(get(varargin{1},'Type'),'axes');
    Ha = varargin{1};
    xtickpos = varargin{2};
    xtickstring = varargin{3};
    if nargin > 3
        textopts = varargin(4:end);
    end
else
    Ha = gca;
    %Hfig = get(Ha,'Parent');
    xtickpos = varargin{1};
    xtickstring = varargin{2};
    if nargin > 2
        textopts = varargin(3:end);
    end
end

set(Ha,'XTick',xtickpos, 'XTickLabel','')
h_olds = findobj(Ha, 'Tag', 'MUXTL');
if ~isempty(h_olds)
    delete(h_olds)
end

%% Make XTickLabels 
NTick = length(xtickpos);
Ybot = max(get(gca,'YLim'));
ht = zeros(NTick,1);
for ii = 1:NTick
    ht(ii) = text('String',xtickstring{ii}, ...
        'Units','data', ...
        'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'center', ...
        'Position',[xtickpos(ii) Ybot], ...
        'Tag','MUXTL');
end
if ~isempty(textopts)
    set(ht,textopts{:})
end

%% shift texts left to fit & squeeze axis if needed

set(Ha,'Units','pixels')
Axpos = get(Ha,'Position');
% set(Hfig,'Units','pixels')
% Figpos = get(Hfig,'Position');

set(ht,'Units','pixels')
TickPos = zeros(NTick,3);
TickExt = zeros(NTick,4);
for ii = 1:NTick
    TickPos(ii,:) = get(ht(ii),'Position');
    tmpext = get(ht(ii),'Extent');
    set(ht(ii),'Position',[TickPos(ii,1) TickPos(ii,2)-tmpext(4)/2-5 TickPos(ii,end)])
    TickExt(ii,:) = get(ht(ii),'Extent');
end

needmove = (Axpos(2) + min(TickExt(:,1)));

if needmove>0;
    Axpos(2) = Axpos(2)+needmove+2;
    Axpos(4) = Axpos(4)-needmove+2;
    set(Ha,'Position',Axpos);
end

set(Ha,'Units','normalized')
set(ht,'Units','normalized')
