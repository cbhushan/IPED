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

/* This function estimates histogram of 1D,2D...ND image using cubic B-spline parzen window estimate.
 * 
 * hist_count = histogram_parzen_variable_size_multithread_double(I, Imin, Imax, nbins, window_size, nthreads);
 *
 * Copyright (c) 2013, C. Bhushan, USC
 */

int mindex2(int x, int y, int sizx) { return y*sizx+x; }

voidthread parzenhistogram(double **Args) {
   double *I1, *Imin, *Imax, *nbins, *scav, *npixels_inv, *ThreadID, *Nthreadsd, *windowSize;
   double *hist_count, *npixels_d;
   int Nthreads, npixels=1;
   
   /* intensity location*/
   int nbins_i;
   double xd, xm, xd_dist, xmd_rel;
   double yd, ym, yd_dist, ymd_rel;
   double fx, fy, segmt_len;
   int i, j, k; /* loop vars*/
   
   I1 = Args[0];
   Imin = Args[1];
   Imax = Args[2];
   nbins = Args[3];
   scav = Args[4];       /* scaling parameter */
   ThreadID = Args[5];
   Nthreadsd = Args[6];
   hist_count = Args[7];
   npixels_d = Args[8];
   windowSize = Args[9];
   
   npixels = (int)npixels_d[0];
   Nthreads = (int)Nthreadsd[0];
   nbins_i = (int)nbins[0];
   segmt_len = 0.25 * windowSize[0];
   
   
   for (i=(int)ThreadID[0]; i<npixels; i=i+Nthreads) {
      
      /* scale image to [0 nbins-1] */
      xd=scav[0]*(I1[i]-Imin[0]);
      
      for (j=ceil(xd-2*segmt_len); j<=floor(xd+2*segmt_len); j++) {
         
         if(j>=0 && j<=(nbins_i-1)) {
            xd_dist = fabs(xd-(double)j);
            
            xmd_rel = xd_dist/segmt_len;
            if(xmd_rel<=1.0) {
               fx = (pow(xmd_rel,3.0)/2.0 - pow(xmd_rel, 2.0) + 2.0/3.0)/segmt_len;
            }
            else {
               fx = (pow(xmd_rel, 2.0) - pow(xmd_rel,3.0)/6.0 - 2.0*xmd_rel + 4.0/3.0)/segmt_len;
            }            
            hist_count[j] += fx;
         }
      }
   }
}

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] ) {
   
   double *I1, *Imin, *Imax, *nbins, *windowSize, *Nthreadsd; /* input */
   double *hist_count; /* outputs*/

   /* For worker threads */
   int Nthreads;
   double ***ThreadArgs, **ThreadArgs1, **ThreadID, *ThreadID1;
   ThreadHANDLE *ThreadList;

   mxArray **histThread;  
   const mwSize *idims; /* Size of input image */
   int odims, nsubs, i, j, npixels=1;
   double npixelsd[1]={1.0};
   double scav[1]={0.0};
   double *histTemp;

   /* Check for proper number of arguments. */
   if(nrhs!=6) {
      mexErrMsgTxt("Six inputs are required.");
   } else if(nlhs!=1) {
      mexErrMsgTxt("Excatly one output is required");
   }
   
   nsubs = mxGetNumberOfDimensions(prhs[0]); /*  Get the number of dimensions */
   idims = mxGetDimensions(prhs[0]); /* Get the sizes of the grid */
   
   /* Calculate total Number of pixels */
   for (i=0; i<nsubs; i++) { 
      npixels=npixels*idims[i]; 
   }
   npixelsd[0]=(double)npixels;
   
   
   /* Assign pointers to each input. */
   I1=(double *)mxGetData(prhs[0]);
   Imin=(double *)mxGetData(prhs[1]);
   Imax=(double *)mxGetData(prhs[2]);
   nbins=(double *)mxGetData(prhs[3]);
   windowSize=(double *)mxGetData(prhs[4]);
   Nthreadsd=(double *)mxGetData(prhs[5]);

   /* Create image matrix for the return arguments*/
   odims=(int)nbins[0]; 
   plhs[0] = mxCreateDoubleMatrix(odims, 1, mxREAL);
   hist_count=(double *)mxGetData(plhs[0]);   /* Assign pointers to output. */
   

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
      histThread[i]=mxCreateDoubleMatrix(odims,1,mxREAL);
      
      /*  Make Thread Structure  */
      ThreadArgs1 = (double **)malloc( 12* sizeof( double * ) );
      ThreadArgs1[0]=I1;
      ThreadArgs1[1]=Imin;
      ThreadArgs1[2]=Imax;
      ThreadArgs1[3]=nbins;
      ThreadArgs1[4]=scav;
      ThreadArgs1[5]=ThreadID[i];
      ThreadArgs1[6]=Nthreadsd;
      ThreadArgs1[7]=(double *)mxGetData(histThread[i]);
      ThreadArgs1[8]=npixelsd;
      ThreadArgs1[9]=windowSize;
      ThreadArgs[i]=ThreadArgs1;
      
      StartThread(ThreadList[i], &parzenhistogram, ThreadArgs[i]) 
   }
   
   for (i=0; i<Nthreads; i++) {
      WaitForThreadFinish(ThreadList[i]);
   }

   for (i=0; i<Nthreads; i++) {
      
      /* add output from each thread */
      histTemp=(double *)mxGetData(histThread[i]);
      for (j=0; j<((int)nbins[0]); j++) {
         hist_count[j] += histTemp[j];
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


