% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 



%%  generate log_lookup.mat 
clc; clear; close all

% for joint pdf
low_lim1 = 1e-9;
up_lim1 = 6e-03;
step1 = low_lim1;

temp1 = low_lim1 : step1 : up_lim1;
log_lookup1 = log(temp1);
range1 = [temp1(1) temp1(end)];

% for marginal 
low_lim2 = 1e-7;
up_lim2 = 0.1;
step2 = low_lim2;

temp2 = low_lim2 : step2 : up_lim2;
log_lookup2 = log(temp2);
range2 = [temp2(1) temp2(end)];

figure; 
subplot(1,2,1); plot(temp1, log_lookup1); xlim([range1(1) range2(2)]); ylim([log_lookup1(1) log_lookup2(end)]); grid on
subplot(1,2,2); plot(temp2, log_lookup2); xlim([range1(1) range2(2)]); ylim([log_lookup1(1) log_lookup2(end)]); grid on

save log_lookup.mat range1 log_lookup1 range2 log_lookup2 step1 step2

%% test accuracy of log lookup table
clc; clear all; close all

tbl = load('log_lookup.mat');
bins = 256;
niter_test = 2000;

% for range1
tbl.md1_rng1 = ones(bins*bins,1)/tbl.range1(1);
tbl.md2_rng1 = ones(2*bins, 1)/tbl.range1(1);

% for range2
tbl.md1_rng2 = ones(bins*bins,1)/tbl.range2(1);
tbl.md2_rng2 = ones(2*bins, 1)/tbl.range2(1);



err = zeros(niter_test, 2);
f_range = 1e-8;

% mode 1
for k = 1:niter_test
      
   % mode = 1, when data is bins x bins (for joint histogram)
   data = rand(bins, bins)*f_range; % 0.125 b/c full range of lookup is [1e-9  0.1]
   d_log = log(data);
   
   % lookup method 
   d_log_lookup = zeros(size(data));
   
   ind_0 = (data<tbl.range1(1));
   ind_1 = (data<tbl.range1(2)) & ~ind_0;
   ind_2 = (data<=tbl.range2(2)) & ~(ind_1 | ind_0);
   ind_3 = (data>tbl.range2(2));
   
   d_log_lookup(ind_0) = log(tbl.range1(1)/2); % constant for lowest quadrant   
   
   ind = uint32(data(ind_1).*tbl.md1_rng1(ind_1));
   d_log_lookup(ind_1) = tbl.log_lookup1(ind);
   
   ind = uint32(data(ind_2).*tbl.md1_rng2(ind_2));
   d_log_lookup(ind_2) = tbl.log_lookup2(ind);
   
   d_log_lookup(ind_3) = log(data(ind_3)); % outside lookup table
   
   err(k,1) = norm(d_log(:)-d_log_lookup(:))/sqrt(numel(data));
   
   
   
   % mode = 2, when data is 2*bins (for concatenated marginal histogram)
   data = rand(2*bins, 1)*f_range; % b/c full range of lookup is [1e-9  0.1]  
   d_log = log(data);
   
   % lookup method 
   d_log_lookup = zeros(size(data));
   
   ind_0 = (data<tbl.range1(1));
   ind_1 = (data<tbl.range1(2)) & ~ind_0;
   ind_2 = (data<=tbl.range2(2)) & ~(ind_1 | ind_0);
   ind_3 = (data>tbl.range2(2));
   
   d_log_lookup(ind_0) = log(tbl.range1(1)/2); % constant for lowest quadrant   
   
   ind = uint32(data(ind_1).*tbl.md2_rng1(ind_1));
   d_log_lookup(ind_1) = tbl.log_lookup1(ind);
   
   ind = uint32(data(ind_2).*tbl.md2_rng2(ind_2));
   d_log_lookup(ind_2) = tbl.log_lookup2(ind);
   
   d_log_lookup(ind_3) = log(data(ind_3)); % outside lookup table

   err(k,2) = norm(d_log(:)-d_log_lookup(:))/sqrt(numel(data));
   
   k
end

figure;
plot(err(:,1));
hold on;
plot(err(:,2), 'r');
title('Norm of error in log estimate')
xlabel('trials')
ylabel('error')
legend('mode = 1', 'mode = 2')
grid on

%    clf;
%    subplot(2,3,1); hist(d(:), 100); title('data')
%    subplot(2,3,2); hist(d_log(:), 100); title('log(data)')
%    subplot(2,3,3); hist(d_log_lookup(:), 100); title('log_lookup(data)')
%
%    subplot(2,3,4); imagesc(abs(d_log-d_log_lookup)); title('abs err'); colorbar
%    subplot(2,3,5); imagesc(d_log); title('log(data)'); colorbar; cax = caxis;
%    subplot(2,3,6); imagesc(d_log_lookup); title('log_lookup(data)'); colorbar; caxis(cax);


%% test speedup  of logFastHandle()
clc; clear all; close all

tbl = load('log_lookup.mat');
bins = 256;
niter_test = 2000;

% for range1
tbl.md1_rng1 = ones(bins*bins,1)/tbl.range1(1);
tbl.md2_rng1 = ones(2*bins, 1)/tbl.range1(1);

% for range2
tbl.md1_rng2 = ones(bins*bins,1)/tbl.range2(1);
tbl.md2_rng2 = ones(2*bins, 1)/tbl.range2(1);



% mode 1
data = rand(bins, bins, niter_test)*1e-17; 

t = tic();
for k = 1:niter_test
   d_log = log(data(:,:,k));
end
t_native = toc(t);

t = tic();
for k = 1:niter_test
   d = data(:,:,k);
   
   d_log_lookup = zeros(size(d));
   
   ind_0 = (d<tbl.range1(1));
   ind_1 = (d<tbl.range1(2)) & ~ind_0;
   ind_2 = (d<=tbl.range2(2)) & ~(ind_1 | ind_0);
   ind_3 = (d>tbl.range2(2));
   
   d_log_lookup(ind_0) = log(tbl.range1(1)/2);
   
   ind = uint32(d(ind_1).*tbl.md1_rng1(ind_1));
   d_log_lookup(ind_1) = tbl.log_lookup1(ind);
   
   ind = uint32(d(ind_2).*tbl.md1_rng2(ind_2));
   d_log_lookup(ind_2) = tbl.log_lookup2(ind);
   
   d_log_lookup(ind_3) = log(d(ind_3));
end
t1 = toc(t);


fprintf('Mode 1 speedup: %f', t_native/t1)

%% real data example
clc; clear; 

dr = '/mnt/work/Projects/BrainReg_Git/data/BDP_rigid_test/single_subject';
bins = 200;
parzen_window_size = 12;
nthreads = 8;

t1 = load_untouch_nii_gz(fullfile(dr, '7151JH.0003.bfc.D_coord.nii.gz'));
epi = load_untouch_nii_gz(fullfile(dr, '7151JH.dwi.00.L0.20.pad.nii.gz'));
msk = load_untouch_nii_gz(fullfile(dr, '7151JH.dwi.00.L0.20.pad.mask.nii.gz'));

t1 = normalize_intensity(t1.img, [0.1 99.99], msk.img>0);
epi = normalize_intensity(epi.img, [0.1 98], msk.img>0);
msk = msk.img>0;

[histVU] = mutual_histogram_parzen_variable_size_multithread_double( double(t1(msk)), double(epi(msk)), ...
                                double(0), double(1), double(bins), ...
                                double(parzen_window_size), double(nthreads) );


 

% for range1
tbl = load('log_lookup.mat');
tbl.md1_rng1 = ones(bins*bins,1)/tbl.range1(1);
tbl.md2_rng1 = ones(2*bins, 1)/tbl.range1(1);

% for range2
tbl.md1_rng2 = ones(bins*bins,1)/tbl.range2(1);
tbl.md2_rng2 = ones(2*bins, 1)/tbl.range2(1);



% mode 1

niter_test = 200;
t = tic();
for k = 1:niter_test
   d_log = log(histVU);
end
t_native = toc(t)

d = histVU;
t = tic();
for k = 1:niter_test
   
   d_log_lookup = zeros(size(d));
   
   ind_0 = (d<tbl.range1(1));
   ind_1 = (d<tbl.range1(2)) & ~ind_0;
   ind_2 = (d<=tbl.range2(2)) & ~(ind_1 | ind_0);
   ind_3 = (d>tbl.range2(2));
   
   d_log_lookup(ind_0) = log(tbl.range1(1)/2);
   
   ind = uint32(d(ind_1).*tbl.md1_rng1(ind_1));
   d_log_lookup(ind_1) = tbl.log_lookup1(ind);
   
   ind = uint32(d(ind_2).*tbl.md1_rng2(ind_2));
   d_log_lookup(ind_2) = tbl.log_lookup2(ind);
   
   d_log_lookup(ind_3) = log(d(ind_3));
end
t1 = toc(t)

fprintf('Mode 1 speedup: %f', t_native/t1)

d = rand(20560)*1e-5;
tic; 
log(d);
toc

d = nan(size(d));
tic; 
log(d);
toc



                             
                             
