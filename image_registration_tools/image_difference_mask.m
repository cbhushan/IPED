% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [t,I]=image_difference_mask(V,U,type,V_mask,U_mask,MaskNum)
% This function gives a registration error and error image I between the two
% images or volumes.
%
% [t,I]=image_difference(I1,I2,type,Mask)
%
% inputs,
%   I1: Input image 1
%   I2: Input image 2
%   type: Type of similarity / error measure
% (optional)
%   Mask: Image/volume which is multiplied with the individual pixel errors
%         before calculation of the te total (mean) similarity error.
%
% if type,
% 'd'  : differences between I1 and I2
% 'sd' : squared differences
% 'mi' : normalized mutual information
% 'mip': normalized mutual information with image split in multiple small regions
% 'gd' : gradient differences
% 'gc' : gradient correlation
% 'cc' : normalized cros correlation
% 'pi' : pattern intensity
% 'ld' : log absolute difference
%
%  Example,
%    I1=im2double(imread('lenag1.png'));
%    I2=im2double(imread('lenag2.png'));
%    [t,I] = image_difference(I1,I2,'sd');
%    disp(t);
%    imshow(I,[])
%
% This function is written by D.Kroon University of Twente (April 2009)

if(exist('type','var')==0), type='sd'; end
%if(exist('Mask','var')==0), Mask=[]; end

if(exist('MaskNum','var')==0), MaskNum=numel(V); end

if(MaskNum==0), t=2; I=2; return; end
if(isempty(V)), t=2; I=2; return; end

% If color image make mask for RGB
% if(~isempty(Mask)&&(size(V,3)==3)&&(size(Mask,3)==1))
% 	Mask_rgb(:,:,1)=Mask; Mask_rgb(:,:,2)=Mask; Mask_rgb(:,:,3)=Mask;
% 	Mask=Mask_rgb;
% end

switch(type)
   case 'd'
      Mask = (V_mask>0) | (U_mask>0);
      if ( nargout > 1 )
         [t,I]=registration_error_differences(V,U,Mask,MaskNum);
      else
         t=registration_error_differences(V,U,Mask,MaskNum);
      end
   case 'sd'
      Mask = (V_mask>0) | (U_mask>0);
      if ( nargout > 1 )
         [t,I]=registration_error_squared_differences(V,U,Mask,MaskNum);
      else
         t=registration_error_squared_differences(V,U,Mask,MaskNum);
      end
      
   case 'mi'
      Mask = (V_mask>0) & (U_mask>0);
      if any(Mask(:))
         if ( nargout > 1 )
            [t,I]=registration_error_mutual_info(V,U,Mask);
         else
            t=registration_error_mutual_info(V,U,Mask);
         end
      else
         t=2; I=2; return;
      end
   case 'mip'
      Mask = (V_mask>0) & (U_mask>0);
      if any(Mask(:))
         if ( nargout > 1 )
            [t,I]=registration_error_local_mutual_info(V,U,Mask);
         else
            t=registration_error_local_mutual_info(V,U,Mask);
         end
      else
         t=2; I=2; return;
      end
   case 'gd'
      Mask = (V_mask>0) | (U_mask>0);
      if ( nargout > 1 )
         [t,I]=registration_error_gradient_difference(V,U,Mask,MaskNum);
         I=2-I;
      else
         t=registration_error_gradient_difference(V,U,Mask,MaskNum);
      end
      t=2-t;
   case 'gc'
      Mask = (V_mask>0) | (U_mask>0);
      if ( nargout > 1 )
         [t,I]=registration_error_gradient_correlation(V,U,Mask,MaskNum);
         I=1-I;
      else
         t=registration_error_gradient_correlation(V,U,Mask,MaskNum);
      end
      t=1-t;
   case 'cc'
      Mask = (V_mask>0) | (U_mask>0);
      if ( nargout > 1 )
         [t,I]=registration_error_normalized_cross_correlation(V,U,Mask,MaskNum);
         I=1-I;
      else
         t=registration_error_normalized_cross_correlation(V,U,Mask,MaskNum);
      end
      t=1-t;
   case 'pi'
      Mask = (V_mask>0) | (U_mask>0);
      if ( nargout > 1 )
         [t,I]=registration_error_pattern_intensity(V,U,Mask,MaskNum);
         I=1-I;
      else
         t=registration_error_pattern_intensity(V,U,Mask,MaskNum);
      end
      t=1-t;
   case 'ld'
      Mask = (V_mask>0) | (U_mask>0);
      if ( nargout > 1 )
         [t,I]=registration_error_log_absolute_distance(V,U,Mask,MaskNum);
      else
         t=registration_error_log_absolute_distance(V,U,Mask,MaskNum);
      end
   otherwise
      error('Unknown error type')
end
if(isnan(t)), warning('imagedifference:NaN','NaN in error image'); t=2; I=2; end


function [t,I]=registration_error_log_absolute_distance(V,U,Mask,MaskNum)
I=log(abs(V-U)+1);
if(~isempty(Mask)), I=I.*Mask;  end
t=sum(I(:))/MaskNum;

function [t,I]=registration_error_normalized_cross_correlation(V,U,Mask,MaskNum)
Vvar=V-mean(V(:)); Uvar=U-mean(U(:));
I=(Vvar.*Uvar)/((sqrt(sum(Vvar(:).^2))*sqrt(sum(Uvar(:).^2)))+eps);
if(~isempty(Mask)), I=I.*Mask;  end
t=sum(I(:))/MaskNum;

function [t,I]=registration_error_gradient_correlation(V,U,Mask,MaskNum)
if(size(U,3)<4)
   [Gx,Gy]=sobel2();
   I=zeros(size(U));
   for i=1:size(U,3)
      if(i==1)
         I(:,:,i)=(1/2)*(registration_error_normalized_cross_correlation(conv2(V(:,:,i),Gx,'same'),conv2(U(:,:,i),Gx,'same'))...
            +registration_error_normalized_cross_correlation(conv2(V(:,:,i),Gy,'same'),conv2(U(:,:,i),Gy,'same')));
      end
   end
else
   [Gx,Gy,Gz]=sobel3();
   I=(1/3)*(registration_error_normalized_cross_correlation(convn(V,Gx,'same'),convn(U,Gx,'same'))...
      +registration_error_normalized_cross_correlation(convn(V,Gy,'same'),convn(U,Gy,'same'))...
      +registration_error_normalized_cross_correlation(convn(V,Gz,'same'),convn(U,Gz,'same')));
end
if(~isempty(Mask)), I=I.*Mask;  end
t=sum(I(:))/MaskNum;

function [Gx,Gy]=sobel2()
Gx=[1 0 -1;2 0 -2;1 0 -1];  Gy=[1 2 1;0 0 0;-1 -2 -1];

function [Gx,Gy,Gz]=sobel3()
Gx=zeros(3,3,3);Gy=zeros(3,3,3); Gz=zeros(3,3,3);
Gx(:,:,1)=[-1 -3 -1;-3 -6 -3;-1 -3 -1]; Gx(:,:,2)=[ 0  0  0; 0  0  0; 0  0  0]; Gx(:,:,3)=[ 1  3  1; 3  6  3; 1  3  1];
Gy(:,:,1)=[ 1  3  1; 0  0  0;-1 -3 -1]; Gy(:,:,2)=[ 3  6  3; 0  0  0;-3 -6 -3]; Gy(:,:,3)=[ 1  3  1; 0  0  0;-1 -3 -1];
Gz(:,:,1)=[-1  0  1;-3  0  3;-1  0  1]; Gz(:,:,2)=[-3  0  3;-6  0  6;-3  0  3]; Gz(:,:,3)=[-1  0  1;-3  0  3;-1  0  1];

function [t,I]=registration_error_squared_differences(V,U,Mask,MaskNum)
if(isempty(Mask)&&(nargout==1))
   if(isa(V,'double'))
      t=squared_difference_double(double(V),double(U));
   else
      t=squared_difference_single(single(V),single(U));
   end
else
   I=(V-U).^2;
   if(~isempty(Mask)), I=I.*Mask;  end
   t=sum(I(:))/(MaskNum);
end

function [t,I]=registration_error_differences(V,U,Mask,MaskNum)
I=(V-U);
if(~isempty(Mask)), I=I.*Mask;  end
t=sum(I(:))/MaskNum;

function [t,I]=registration_error_pattern_intensity(V,U,Mask,MaskNum)
Idiff=V./(mean(V(:))+1e-5)-U./(mean(U(:))+1e-5);
o=0.3; %determines if grey-value varion is a structure (must be laster than noise)
r=5; numr=0; listr=[];
if(size(U,3)<4)
   if(size(U,3)>1), U=rgb2gray(U); V=rgb2gray(V); Idiff=V./(mean(V(:))+1e-5)-U./(mean(U(:))+1e-5); end
   
   for u=-r:r
      for v=-r:r
         if((u^2+v^2)<=r^2), numr=numr+1; listr(numr,:)=[u v]; end
      end
   end
   for u=-r:r
      for v=-r:r
         if((u^2+v^2)<=r^2), numr=numr+1; listr(numr,:)=[u v]; end
      end
   end
   SP1=zeros(size(U)-2*r,class(Idiff));
   for i=1:size(listr,1),
      u=listr(i,1); v=listr(i,2);
      SP1=SP1+o^2./(o^2+(Idiff(1+r:end-r,1+r:end-r)-Idiff(1+r+u:end-r+u,1+r+v:end-r+v)).^2);
   end
else
   for u=-r:r,
      for v=-r:r
         for w=-r:r
            if((u^2+v^2+w^2)<=r^2), numr=numr+1; listr(numr,:)=[u v w]; end
         end
      end
   end
   SP1=zeros(size(U)-2*r,class(Idiff));
   for i=1:size(listr,1),
      u=listr(i,1); v=listr(i,2); w=listr(i,3);
      SP1=SP1+o^2./(o^2+(Idiff(1+r:end-r,1+r:end-r,1+r:end-r)-Idiff(1+r+u:end-r+u,1+r+v:end-r+v,1+r+w:end-r+w)).^2);
   end
end
I=SP1./(size(listr,1)+eps);
if(~isempty(Mask)), I=I.*Mask;  end
t=sum(I(:))/MaskNum;


function [t,I]=registration_error_gradient_difference(V,U,Mask,MaskNum)
if(size(U,3)<4)
   Sgdiff=zeros(size(U));
   for i=1:size(U,3)
      a=mean(V(:))/mean(U(:));
      [Gx,Gy]=sobel2();
      Idiffv=conv2(V(:,:,i),Gx,'same')-a*conv2(U(:,:,i),Gx,'same');
      Idiffh=conv2(V(:,:,i),Gy,'same')-a*conv2(U(:,:,i),Gy,'same');
      Av=var(Idiffv(:));
      Ah=var(Idiffh(:));
      Sgdiff(:,:,i)=(Av./(Av+Idiffv.^2+eps))+(Ah./(Ah+Idiffh.^2+eps));
   end
else
   a=mean(V(:))/mean(U(:));
   [Gx,Gy,Gz]=sobel3();
   Idiffv=convn(V,Gx,'same')-a*convn(U,Gx,'same');
   Idiffh=convn(V,Gy,'same')-a*convn(U,Gy,'same');
   Idiffz=convn(V,Gz,'same')-a*convn(U,Gz,'same');
   Av=var(Idiffv(:));
   Ah=var(Idiffh(:));
   Az=var(Idiffz(:));
   Sgdiff=(Av./(Av+Idiffv.^2+eps))+(Ah./(Ah+Idiffh.^2+eps))+(Az./(Az+Idiffz.^2+eps));
end
I=Sgdiff;
if(~isempty(Mask)), I=I.*Mask;  end
t=sum(I(:))/MaskNum;


function [t,I]=registration_error_local_mutual_info(V,U,Mask)
% Split the image in multiple regions to allow a more local mutual information measure.

sub_size = round(size(V)/2);
lim = [1 1 1; sub_size; round(sub_size*0.5); round(sub_size*0.5)+sub_size-1;  size(V)-sub_size+1; size(V)];

Vpart = zeros([sub_size 27]);
Upart = zeros([sub_size 27]);
Maskpart = zeros([sub_size 27]);

counter = 1;
for x=1:2:6
   x_lim = lim(x,1):lim(x+1,1);
   for y=1:2:6
      y_lim = lim(y,2):lim(y+1,2);
      for z=1:2:6
         z_lim = lim(z,3):lim(z+1,3);
         Vpart(:,:,:,counter) = V(x_lim, y_lim, z_lim);
         Upart(:,:,:,counter) = U(x_lim, y_lim, z_lim);
         Maskpart(:,:,:,counter) = Mask(x_lim, y_lim, z_lim);
         counter = counter+1;
      end
   end
end

t = 0;
for k = 1:27
   if any(Maskpart(:,:,:,k)>0)
      tpart = registration_error_mutual_info(Vpart(:,:,:,k),Upart(:,:,:,k),Maskpart(:,:,:,k)>0);
      t = t+tpart;
   end
end
I=[];

function [t,I]=registration_error_mutual_info(V,U,Mask)
% registration error based on normalized mutual information. Estimates
% joint pdf  using cubic b-spline parzen window.
%
%      err = (H(V) + H(U)) / H(V,U)
% 

global bins nthreads;%=128;

% Remove unmasked pixels
Vm = V(Mask);
Um = U(Mask);

%range = [min([Vm(:);Um(:)])  max([Vm(:);Um(:)])];
range = [0 1];

[histVU] = mutual_histogram_parzen_multithread_double(double(Vm), double(Um), double(range(1)), double(range(2)), double(bins), double(nthreads));

histV = double(sum(histVU, 1));
histU = double(sum(histVU, 2));
histVU = double(histVU);

% Calculate fast log (lookup table)
global log_lookup1 scale1 range1 ;
global log_lookup2 scale2 range2 ;

log_histVU = zeros(size(histVU));
ind_0 = (histVU<range1(1));
log_histVU(ind_0) = log(range1(1));
ind_1 = (histVU<=range1(2)) & ~ind_0;
ind_2 = (histVU>range1(2));
ind = histVU(ind_1).*scale1(ind_1);
log_histVU(ind_1) = log_lookup1(uint32(ind));
log_histVU(ind_2) = log(histVU(ind_2)); % outside lookup table

histV_histU = [histV(:); histU(:)];
log_histV_histU = zeros(size(histV_histU));
ind_0 = (histV_histU<range2(1));
log_histV_histU(ind_0) = log(range2(1));
ind_1 = (histV_histU<=range2(2)) & ~ind_0;
ind_2 = (histV_histU>range2(2));
ind = histV_histU(ind_1).*scale2(ind_1);
log_histV_histU(ind_1) = log_lookup2(uint32(ind));
log_histV_histU(ind_2) = log(histV_histU(ind_2)); % outside lookup table

log_histV = log_histV_histU(1:bins);
log_histU = log_histV_histU(bins+1:end);

p1log = histV(:) .* log_histV(:);
p2log = histU(:) .* log_histU(:);
p12log = histVU .* log_histVU;

% Calculate amount of Information
MI_num = sum(p1log) + sum(p2log);
MI_den = sum(p12log(:));

% Studholme, Normalized mutual information
if(MI_den==0)
   MI_den = eps;   
   error('check this case.') 
end

t = -1*(MI_num/MI_den);
I = [];
