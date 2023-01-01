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
 * Modified version of mutual_histogram_parzen_double.c to include multi-threading 
 * and variable size parzen window.
 *
 * hist12 = mutual_histogram_parzen_variable_size_multithread_double(I1, I2, Imin, Imax, nbins, window_size, nthreads);
 * hist1 = sum(hist12, 1);
 * hist2 = sum(hist12, 2);
 *
 * Original license below: 
 * http://www.mathworks.com/matlabcentral/fileexchange/21451-multimodality-non-rigid-demon-algorithm-image-registration
 * 
 * Copyright (c) 2009, Dirk-Jan Kroon
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in
 *   the documentation and/or other materials provided with the distribution
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

int mindex2(int x, int y, int sizx) { return y*sizx+x; }

voidthread parzenhistogram(double **Args) {
   double *I1, *I2, *Imin, *Imax, *nbins, *scav, *npixels_inv, *ThreadID, *Nthreadsd, *windowSize;
   double *hist12, *npixels_d;
   int Nthreads, npixels=1;

   /* intensity location*/
   int nbins_i;
   double xd, xm, xd_dist, xmd_rel;
   double yd, ym, yd_dist, ymd_rel;
   double fx, fy, segmt_len;
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
   windowSize = Args[11];
   
   npixels = (int)npixels_d[0];
   Nthreads = (int)Nthreadsd[0];
   nbins_i=(int) nbins[0];
   segmt_len = 0.25 * windowSize[0];
   
   
   for (i=(int)ThreadID[0]; i<npixels; i=i+Nthreads) {
      
      /* scale image to [0 nbins-1] */
      xd=scav[0]*(I1[i]-Imin[0]); 
      yd=scav[0]*(I2[i]-Imin[0]); 
      
      for (j=ceil(xd-2*segmt_len); j<=floor(xd+2*segmt_len); j++) {
         for (k=ceil(yd-2*segmt_len); k<=floor(yd+2*segmt_len); k++) {

            if(j>=0 && j<=(nbins_i-1) && k>=0 && k<=(nbins_i-1)) {
               xd_dist = fabs(xd-(double)j);
               yd_dist = fabs(yd-(double)k);
               
               xmd_rel = xd_dist/segmt_len;
               if(xmd_rel<=1.0) {
                  fx = (pow(xmd_rel,3.0)/2.0 - pow(xmd_rel, 2.0) + 2.0/3.0)/segmt_len;
               }
               else {
                  fx = (pow(xmd_rel, 2.0) - pow(xmd_rel,3.0)/6.0 - 2.0*xmd_rel + 4.0/3.0)/segmt_len;
               }
               
               ymd_rel = yd_dist/segmt_len;
               if(ymd_rel<=1.0) {
                  fy = (pow(ymd_rel,3.0)/2.0 - pow(ymd_rel, 2.0) + 2.0/3.0)/segmt_len;
               }
               else {
                  fy = (pow(ymd_rel, 2.0) - pow(ymd_rel,3.0)/6.0 - 2.0*ymd_rel + 4.0/3.0)/segmt_len;
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
   
   double *I1, *I2, *Imin, *Imax, *nbins, *windowSize, *Nthreadsd; /* input */
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
   if(nrhs!=7) {
      mexErrMsgTxt("seven inputs are required.");
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
   windowSize=(double *)mxGetData(prhs[5]);
   Nthreadsd=(double *)mxGetData(prhs[6]);

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
      ThreadArgs1 = (double **)malloc( 12* sizeof( double * ) );
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
      ThreadArgs1[11]=windowSize;
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


