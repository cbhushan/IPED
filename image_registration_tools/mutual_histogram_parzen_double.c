/*
 * 
 * IPED - Improved B0-distortion correction in diffusion MRI
 * Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
 * 
 * This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
 * Please see https://github.com/cbhushan/IPED for details.
 * 
 * SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
 * 
 */


#include "mex.h"
#include "math.h"

/* This function estimates 2D joint histogram of 1D,2D...ND images
 * using cubic B-spline parzen window estimate.
 *
 * hist12 = mutual_histogram_double(I1,I2,Imin,Imax,nbins);
 * hist1 = sum(hist12, 1);
 * hist2 = sum(hist12, 2);
 */

int mindex2(int x, int y, int sizx) { return y*sizx+x; }

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
   /*   I1 and I2 are the input images */
   /*  hist12 joint histogram */
   double *I1, *I2, *Imin, *Imax, *nbins, *hist12;
   
   /* Size of input image */
   const mwSize *idims;
   
   /*  Size of output */
   int odims2[2]={0,0};
   int odims1[1]={0};
   /*  Dimensions */
   int nsubs;
   int npixels=1;
   double npixelsi=1;
   
   /* intensity location*/
   int sizex, sizexd;
   double xd, xm, xp, xmd, xpd;
   double yd, ym, yp, ymd, ypd;
   double fx, fy;
   int xmi, xpi, ymi, ypi;
   
   /* loop vars*/
   int i, j, k;
   
   /*  vars*/
   double minv;
   double scav;
   
   /* Check for proper number of arguments. */
   if(nrhs!=5) {
      mexErrMsgTxt("five inputs are required.");
   } else if(nlhs!=1) {
      mexErrMsgTxt("One outputs are required");
   }
   
   /*  Get the number of dimensions */
   nsubs = mxGetNumberOfDimensions(prhs[0]);
   /* Get the sizes of the grid */
   idims = mxGetDimensions(prhs[0]);
   for (i=0; i<nsubs; i++) { npixels=npixels*idims[i]; }
   npixelsi=1.0/((double)npixels);
   
   /* Assign pointers to each input. */
   I1=(double *)mxGetData(prhs[0]);
   I2=(double *)mxGetData(prhs[1]);
   Imin=(double *)mxGetData(prhs[2]);
   Imax=(double *)mxGetData(prhs[3]);
   nbins=(double *)mxGetData(prhs[4]);
   
   /*  Create image matrix for the return arguments*/
   odims2[0]=(int) nbins[0]; odims2[1]=(int)nbins[0];
   plhs[0] = mxCreateNumericArray(2, odims2, mxDOUBLE_CLASS, mxREAL);
   
   /* Assign pointers to each output. */
   hist12=(double *)mxGetData(plhs[0]);

   /* min value */
   minv=Imin[0];
   /* scale value */
   scav=(nbins[0]-1.0)/(Imax[0]-Imin[0]);
   sizex=(int) nbins[0];
   sizexd=sizex-1;
   
   for (i=0; i<npixels; i++) {
      xd=(double)scav*(I1[i]-minv);
      xm=(double)floor(xd);
      xmi=(int)xm;
      
      yd=(double)scav*(I2[i]-minv);
      ym=(double)floor(yd);
      ymi=(int)ym;
      
      /* putting things outside range on the edge of the bins OR normalize image before input */
      /*if(xd<-1){ xd=-1; xmi=-1; } else if(xd>sizexd+1) { xd=sizexd+1; xmi=sizexd+1; }
      /*if(yd<-1){ yd=-1; ymi=-1; } else if(yd>sizexd+1) { yd=sizexd+1; ymi=sizexd+1; }
      */
      
      for (j=xmi-1; j<xmi+3; j++) {
         for (k=ymi-1; k<ymi+3; k++) {

            if(j>=0 && j<=sizexd && k>=0 && k<=sizexd) {
               xmd = fabs(xd-(double)j);
               ymd = fabs(yd-(double)k);
 
               
               if(xmd<=1.0) {
                  fx = (pow(xmd,3)/2.0 - pow(xmd, 2) + 2.0/3.0);
               }
               else {
                  fx = (pow(xmd, 2) - pow(xmd,3)/6.0 - 2.0*xmd + 4.0/3.0);
               }
               
               if(ymd<=1.0) {
                  fy = (pow(ymd,3)/2.0 - pow(ymd, 2) + 2.0/3.0);
               }
               else {
                  fy = (pow(ymd, 2) - pow(ymd,3)/6.0 - 2.0*ymd + 4.0/3.0);
               }

               hist12[j*sizex + k] += (fx*fy*npixelsi);
            }
         }
      }
   }
}


