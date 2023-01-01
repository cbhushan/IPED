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
#include "multiple_os_thread.h"

/* This function estimates 2D joint histogram of 1D,2D...ND images
 * using cubic B-spline parzen window estimate.
 * 
 * Modified version of mutual_histogram_parzen_variable_size_multithread_double.c to allow specifying different
 * parzen width along different dimensions.
 *
 * hist12 = mutual_histogram_parzen_two_variable_size_multithread_double(I1, I2, Imin, Imax, nbins, window_size1, window_size2, nthreads);
 * hist1 = sum(hist12, 1);
 * hist2 = sum(hist12, 2);
 */

int mindex2(int x, int y, int sizx) { return y*sizx+x; }

voidthread parzenhistogram(double **Args) {
   double *I1, *I2, *Imin, *Imax, *nbins, *scav, *npixels_inv, *ThreadID, *Nthreadsd, *windowSize1, *windowSize2;
   double *hist12, *npixels_d;
   int Nthreads, npixels=1;

   /* intensity location*/
   int nbins_i;
   double xd, xm, xd_dist, xmd_rel;
   double yd, ym, yd_dist, ymd_rel;
   double fx, fy, segmt_len1, segmt_len2;
   int i, j, k; /* loop vars*/
   
   I1 = Args[0];
   I2 = Args[1];
   Imin = Args[2];
   Imax = Args[3];
   nbins = Args[4];
   scav = Args[5];       /* scaling parameter */
   npixels_inv = Args[6];
   ThreadID = Args[7];
   Nthreadsd = Args[8];
   hist12 = Args[9];
   npixels_d = Args[10];
   windowSize1 = Args[11];
   windowSize2 = Args[12];
   
   npixels = (int)npixels_d[0];
   Nthreads = (int)Nthreadsd[0];
   nbins_i=(int) nbins[0];
   segmt_len1 = 0.25 * windowSize1[0];
   segmt_len2 = 0.25 * windowSize2[0];
   
   for (i=(int)ThreadID[0]; i<npixels; i=i+Nthreads) {
      
      /* scale image to [0 nbins-1] */
      xd=scav[0]*(I1[i]-Imin[0]); 
      yd=scav[0]*(I2[i]-Imin[0]); 
      
      for (j=ceil(xd-2*segmt_len1); j<=floor(xd+2*segmt_len1); j++) {
         for (k=ceil(yd-2*segmt_len2); k<=floor(yd+2*segmt_len2); k++) {

            if(j>=0 && j<=(nbins_i-1) && k>=0 && k<=(nbins_i-1)) {
               xd_dist = fabs(xd-(double)j);
               yd_dist = fabs(yd-(double)k);
               
               xmd_rel = xd_dist/segmt_len1;
               if(xmd_rel<=1.0) {
                  fx = (pow(xmd_rel,3.0)/2.0 - pow(xmd_rel, 2.0) + 2.0/3.0)/segmt_len1;
               }
               else {
                  fx = (pow(xmd_rel, 2.0) - pow(xmd_rel,3.0)/6.0 - 2.0*xmd_rel + 4.0/3.0)/segmt_len1;
               }
               
               ymd_rel = yd_dist/segmt_len2;
               if(ymd_rel<=1.0) {
                  fy = (pow(ymd_rel,3.0)/2.0 - pow(ymd_rel, 2.0) + 2.0/3.0)/segmt_len2;
               }
               else {
                  fy = (pow(ymd_rel, 2.0) - pow(ymd_rel,3.0)/6.0 - 2.0*ymd_rel + 4.0/3.0)/segmt_len2;
               }

               hist12[j*nbins_i + k] += (fx*fy*npixels_inv[0]);
            }
         }
      }
   }
}

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] ) {
   
   double *I1, *I2, *Imin, *Imax, *nbins, *windowSize1, *windowSize2, *Nthreadsd; /* input */
   double *hist12; /* outputs*/

   /* For worker threads */
   int Nthreads;
   double ***ThreadArgs, **ThreadArgs1, **ThreadID, *ThreadID1;
   ThreadHANDLE *ThreadList;

   mxArray **histThread;  
   const mwSize *idims; /* Size of input image */
   int odims, nsubs, i, j, npixels=1;
   double npixels_inv[1]={1.0};
   double npixelsd[1]={1.0};
   double scav[1]={0.0};
   double *histTemp;

   /* Check for proper number of arguments. */
   if(nrhs!=8) {
      mexErrMsgTxt("Eight inputs are required.");
   } else if(nlhs!=1) {
      mexErrMsgTxt("Excatly one output is required");
   }
   
   nsubs = mxGetNumberOfDimensions(prhs[0]); /*  Get the number of dimensions */
   idims = mxGetDimensions(prhs[0]); /* Get the sizes of the grid */
   
   /* Calculate total Number of pixels */
   for (i=0; i<nsubs; i++) { 
      npixels=npixels*idims[i]; 
   }
   
   npixels_inv[0]=1.0/((double)npixels);
   npixelsd[0]=(double)npixels;
   
   
   /* Assign pointers to each input. */
   I1=(double *)mxGetData(prhs[0]);
   I2=(double *)mxGetData(prhs[1]);
   Imin=(double *)mxGetData(prhs[2]);
   Imax=(double *)mxGetData(prhs[3]);
   nbins=(double *)mxGetData(prhs[4]);
   windowSize1=(double *)mxGetData(prhs[5]);
   windowSize2=(double *)mxGetData(prhs[6]);
   Nthreadsd=(double *)mxGetData(prhs[7]);

   /*  Create image matrix for the return arguments*/
   odims=(int)nbins[0]; 
   plhs[0] = mxCreateDoubleMatrix(odims, odims, mxREAL);
   hist12=(double *)mxGetData(plhs[0]);   /* Assign pointers to output. */
   

   /* Reserve room for handles of threads in ThreadList  */
   Nthreads=(int)Nthreadsd[0];
   ThreadList = (ThreadHANDLE*)malloc(Nthreads* sizeof( ThreadHANDLE ));
   ThreadID = (double **)malloc( Nthreads* sizeof(double *) );
   ThreadArgs = (double ***)malloc( Nthreads* sizeof(double **) );

   histThread = (mxArray **)malloc(Nthreads* sizeof(mxArray *));
   
   /* scaling parameter */
   scav[0]=(nbins[0]-1.0)/(Imax[0]-Imin[0]);   

   /* Threads */
   for (i=0; i<Nthreads; i++) {
      
      /*  Make Thread ID  */
      ThreadID1= (double *)malloc( 1* sizeof(double) );
      ThreadID1[0]=(double)i;
      ThreadID[i]=ThreadID1;

      /* output of each thread */
      histThread[i]=mxCreateDoubleMatrix(odims,odims,mxREAL);
      
      /*  Make Thread Structure  */
      ThreadArgs1 = (double **)malloc( 13* sizeof( double * ) );
      ThreadArgs1[0]=I1;
      ThreadArgs1[1]=I2;
      ThreadArgs1[2]=Imin;
      ThreadArgs1[3]=Imax;
      ThreadArgs1[4]=nbins;
      ThreadArgs1[5]=scav;
      ThreadArgs1[6]=npixels_inv;
      ThreadArgs1[7]=ThreadID[i];
      ThreadArgs1[8]=Nthreadsd;
      ThreadArgs1[9]=(double *)mxGetData(histThread[i]);
      ThreadArgs1[10]=npixelsd;
      ThreadArgs1[11]=windowSize1;
      ThreadArgs1[12]=windowSize2;
      ThreadArgs[i]=ThreadArgs1;
      
      StartThread(ThreadList[i], &parzenhistogram, ThreadArgs[i]) 
   }
   
   for (i=0; i<Nthreads; i++) {
      WaitForThreadFinish(ThreadList[i]);
   }

   for (i=0; i<Nthreads; i++) {
      
      /* add output from each thread */
      histTemp=(double *)mxGetData(histThread[i]);
      for (j=0; j<((int)nbins[0]*(int)nbins[0]); j++) {
         hist12[j] += histTemp[j];
      }
      
      mxDestroyArray(histThread[i]);
      free(ThreadArgs[i]);
      free(ThreadID[i]);
   }
   
   free(histThread);
   free(ThreadArgs);
   free(ThreadID );
   free(ThreadList);
}


