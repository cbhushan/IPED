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
 * Multi-threaded version of mutual_histogram_parzen_double.c
 *
 * hist12 = mutual_histogram_parzen_multithread_double(I1,I2,Imin,Imax,nbins,nthreads);
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
   double *I1, *I2, *Imin, *Imax, *nbins, *scav, *npixelsi, *ThreadID, *Nthreadsd;
   double *hist12, *npixelsd;
   int Nthreads, npixels=1;

   /* intensity location*/
   int sizex, sizexd;
   double xd, xm, xp, xmd, xpd;
   double yd, ym, yp, ymd, ypd;
   double fx, fy;
   int xmi, xpi, ymi, ypi;
   
   /* loop vars*/
   int i, j, k;
   double minv;
   
   I1 = Args[0];
   I2 = Args[1];
   Imin = Args[2];
   Imax = Args[3];
   nbins = Args[4];
   scav = Args[5];
   npixelsi = Args[6];
   ThreadID = Args[7];
   Nthreadsd = Args[8];
   hist12 = Args[9];
   npixelsd = Args[10];
   
   npixels = (int)npixelsd[0];
   Nthreads = (int)Nthreadsd[0];

   sizex=(int) nbins[0];
   sizexd=sizex-1;
   /* min value */
   minv=Imin[0];
   
   for (i=(int)ThreadID[0]; i<npixels; i=i+Nthreads) {
      xd=scav[0]*(I1[i]-Imin[0]);
      xm=(double)floor(xd);
      xmi=(int)xm;
      
      yd=scav[0]*(I2[i]-Imin[0]);
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
                  fx = (pow(xmd,3.0)/2.0 - pow(xmd, 2.0) + 2.0/3.0);
               }
               else {
                  fx = (pow(xmd, 2.0) - pow(xmd,3.0)/6.0 - 2.0*xmd + 4.0/3.0);
               }
               
               if(ymd<=1.0) {
                  fy = (pow(ymd,3.0)/2.0 - pow(ymd, 2.0) + 2.0/3.0);
               }
               else {
                  fy = (pow(ymd, 2.0) - pow(ymd,3.0)/6.0 - 2.0*ymd + 4.0/3.0);
               }

               hist12[j*sizex + k] += (fx*fy*npixelsi[0]);
            }
         }
      }
   }
}

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] ) {
   
   double *I1, *I2, *Imin, *Imax, *nbins, *Nthreadsd; /* input */
   double *hist12; /* outputs*/

   /* For worker threads */
   int Nthreads;
   double ***ThreadArgs, **ThreadArgs1, **ThreadID, *ThreadID1;
   ThreadHANDLE *ThreadList;

   mxArray **histThread;  
   const mwSize *idims; /* Size of input image */
   int odims, nsubs, i, j, npixels=1;
   double npixelsi[1]={1.0};
   double npixelsd[1]={1.0};
   double scav[1]={0.0};
   double *histTemp;

   /* Check for proper number of arguments. */
   if(nrhs!=6) {
      mexErrMsgTxt("six inputs are required.");
   } else if(nlhs!=1) {
      mexErrMsgTxt("One output is required");
   }
   
   nsubs = mxGetNumberOfDimensions(prhs[0]); /*  Get the number of dimensions */
   idims = mxGetDimensions(prhs[0]); /* Get the sizes of the grid */
   
   for (i=0; i<nsubs; i++) { npixels=npixels*idims[i]; }
   npixelsi[0]=1.0/((double)npixels);
   npixelsd[0]=(double)npixels;
   
   
   /* Assign pointers to each input. */
   I1=(double *)mxGetData(prhs[0]);
   I2=(double *)mxGetData(prhs[1]);
   Imin=(double *)mxGetData(prhs[2]);
   Imax=(double *)mxGetData(prhs[3]);
   nbins=(double *)mxGetData(prhs[4]);
   Nthreadsd=(double *)mxGetData(prhs[5]);
   
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

   scav[0]=(nbins[0]-1.0)/(Imax[0]-Imin[0]);   /* scaling parameter */

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
      ThreadArgs1[6]=npixelsi;
      ThreadArgs1[7]=ThreadID[i];
      ThreadArgs1[8]=Nthreadsd;
      ThreadArgs1[9]=(double *)mxGetData(histThread[i]);
      ThreadArgs1[10]=npixelsd;
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


