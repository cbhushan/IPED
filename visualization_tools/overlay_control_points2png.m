% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function overlay_control_points2png(epiImgR_plot, O_trans, spacing_plot, fname)
% Plot controls points slice wise over image

x = O_trans(:,:,:,1)+1;
y = O_trans(:,:,:,2)+1;

x_max = max(x(:))+4;
x_min = min(x(:))-4;
y_max = max(y(:))+4;
y_min = min(y(:))-4;

h = figure('Units','pixels','Position',[10 10 1000 1000]);
for n = 1:size(O_trans, 3)
   imagesc(epiImgR_plot(:,:,(n-1)*spacing_plot(3)+1), [0 1]);
   colormap(gray)
   hold on
   
   x = O_trans(:,:,n,1)+1;
   y = O_trans(:,:,n,2)+1;

   ylim([x_min x_max]);
   xlim([y_min y_max]);
   set(gca,'Color',[0 0 0.7]);
   axis equal 
   
   plot(y(:),x(:),'.', 'MarkerEdgeColor','r', 'MarkerSize',8);
   set(gcf, 'InvertHardCopy', 'off');
   saveas(gcf, [fname '.' num2str(n) '.png'])
   clf
   hold off
end
close(h)

end
