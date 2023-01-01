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
#include "float.h"

/* 
 * This function makes a 2D joint histogram of 1D,2D...ND images
 * and also calculates the seperate histograms of both images and returns log of the histograms.
 *
 * This function is written to speed-up log calculation by computing log in C. Does NOT help!
 * Use lookup tables instead.
 *
 * [hist12, hist1, hist2]=mutual_histogram_double(I1,I2,Imin,Imax,nbins);
 *
 */

int mindex2(int x, int y, int sizx) { return y*sizx+x; }

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    /*   I1 and I2 are the input images */
    /*  hist12 joint histogram */
    /*  hist1 and histogram of I1 and hist2 of I2 */
    double *I1, *I2, *Imin, *Imax, *nbins, *hist12, *hist12log, *hist1, *hist2, *hist1log, *hist2log;
    
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
    int xmi, xpi, ymi, ypi;
	
    /* loop vars*/
    int i;
    
    /*  vars*/
    double minv;
    double scav;
    
    /* Check for proper number of arguments. */
    if(nrhs!=5) {
       mexErrMsgTxt("five inputs are required.");
    } else if(nlhs!=6) {
       mexErrMsgTxt("Three outputs are required");
    }
  
    /*  Get the number of dimensions */
    nsubs = mxGetNumberOfDimensions(prhs[0]);
    /* Get the sizes of the grid */
    idims = mxGetDimensions(prhs[0]);   
    for (i=0; i<nsubs; i++) { npixels=npixels*idims[i]; }
    npixelsi=1/((double)npixels);

    /* Assign pointers to each input. */
    I1=(double *)mxGetData(prhs[0]);
    I2=(double *)mxGetData(prhs[1]);
    Imin=(double *)mxGetData(prhs[2]);
    Imax=(double *)mxGetData(prhs[3]);
    nbins=(double *)mxGetData(prhs[4]);
    
    /*  Create image matrix for the return arguments*/
    odims2[0]=(int) nbins[0]; odims2[1]=(int)nbins[0];  
    plhs[0] = mxCreateNumericArray(2, odims2, mxDOUBLE_CLASS, mxREAL);

    odims1[0]=(int) nbins[0]; 
    plhs[1] = mxCreateNumericArray(1, odims1, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(1, odims1, mxDOUBLE_CLASS, mxREAL);
    plhs[3] = mxCreateNumericArray(2, odims2, mxDOUBLE_CLASS, mxREAL);
    plhs[4] = mxCreateNumericArray(1, odims1, mxDOUBLE_CLASS, mxREAL);
    plhs[5] = mxCreateNumericArray(1, odims1, mxDOUBLE_CLASS, mxREAL);
    
    /* Assign pointers to each output. */
    hist12=(double *)mxGetData(plhs[0]);
    hist1=(double *)mxGetData(plhs[1]);
    hist2=(double *)mxGetData(plhs[2]);
    hist12log=(double *)mxGetData(plhs[3]);
    hist1log=(double *)mxGetData(plhs[4]);
    hist2log=(double *)mxGetData(plhs[5]);
    
    /* min value */
    minv=Imin[0];
    /* scale value */
    scav=nbins[0]/(Imax[0]-Imin[0]);
    sizex=(int) nbins[0];
    sizexd=sizex-1;

    for (i=0; i<npixels; i++)
    {
        xd=(double)scav*(I1[i]-minv);
        xm=(double)floor(xd); xp=xm+1;
        xmd=xp-xd; xpd=xd-xm;
                
        yd=(double)scav*(I2[i]-minv);
        ym=(double)floor(yd); yp=ym+1;
        ymd=yp-yd; ypd=yd-ym;

        xmi=(int)xm; xpi=(int)xp;
		ymi=(int)ym; ypi=(int)yp;
		
        /* Make sum of all values in histogram 1 and histrogram 2 equal to 1*/
         
        xmd*=npixelsi; ymd*=npixelsi; xpd*=npixelsi;  ypd*=npixelsi;
                        

        if(xmi<0){ xmi=0; } else if(xmi>sizexd) { xmi=sizexd; }
        if(xpi<0){ xpi=0; } else if(xpi>sizexd) { xpi=sizexd; }
        if(ymi<0){ ymi=0; } else if(ymi>sizexd) { ymi=sizexd; }
        if(ypi<0){ ypi=0; } else if(ypi>sizexd) { ypi=sizexd; }

        hist12[xmi+ymi*sizex]+=xmd*ymd;
        hist12[xpi+ymi*sizex]+=xpd*ymd;
        hist12[xmi+ypi*sizex]+=xmd*ypd;
        hist12[xpi+ypi*sizex]+=xpd*ypd;

        hist1[xmi]=hist1[xmi]+xmd; hist1[xpi]=hist1[xpi]+xpd;
        hist2[ymi]=hist2[ymi]+ymd; hist2[ypi]=hist2[ypi]+ypd;
    }
    
    /*printf("DBL_MIN = %f\n", DBL_MIN);*/
    for (i=0; i<((int)nbins[0]*(int)nbins[0]); i++)
    {
       hist12log[i] = log(hist12[i]+DBL_MIN);
    }
    
    for (i=0; i<(int)nbins[0]; i++)
    {
       hist1log[i] = log(hist1[i]+DBL_MIN);
       hist2log[i] = log(hist2[i]+DBL_MIN);
    }
    
}
        

