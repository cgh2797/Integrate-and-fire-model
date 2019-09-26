/*=================================================================
 *
 * rundynam_gif_mex.C	
 *   .MEX code for running dynamics of generalized integrate-and-firemodel 
 *
 * The calling syntax is:
 *
 *   [Sp,Vmem,Ispkhist] = rundynam_gif_mex(Iinj,ih,vleak,vthr,vreset, ...
 *                               sig,dcay1,dcay2,nstep,rndseed)
 *     
 *     Inputs:  Iinj = injected current (column vector), 
 *              ih = h after-current (sampled on same time lattice);
 *              vleak = reversal potential
 *              dcay1 = decay of V during one time bin
 *              dcay2 = integrated decay (i.e. how much injected
 *                      DC current accumulates during a time bin) 
 *              nstep = # times to subsample input current (integer, =1/dtsim)
 *              rndseed = seed for random number generator
 *
 *     Outputs: Sp = binary vector containing spike train
 *              Vmem = voltage trace
 *              Ispkhist = injected current from spike histories
 *
 *    Dependency: gwnoise.cpp
 *
 *  JP  6/19/07
 *=================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "gwnoise.cpp"

size_t sf=sizeof(float), sd=sizeof(double), si=sizeof(long);

/* seed  */
long idum;

/* subroutines  */
float gasdev(long *idum);

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *Iinj, *ih, vleak, sig, dcay1 , dcay2, vthr, vr, *ptr;
  double *Sp, *V, *I, Iprev, vprev, dtsim, II;
  int slen, rlen, nh, ii, jj, kk, j, nstep;
  float nse;
  
  /*  Read inputs into appropriate local variables */
  slen = mxGetM(prhs[0]);
  Iinj = mxGetPr(prhs[0]);  /* injected current */
  nh = mxGetM(prhs[1]);
  ih = mxGetPr(prhs[1]);    /* h current (post-spike) */
  
  ptr = mxGetPr(prhs[2]);
  vleak = ptr[0];
  ptr = mxGetPr(prhs[3]);
  vthr = ptr[0];
  ptr = mxGetPr(prhs[4]);
  vr = ptr[0];
  ptr = mxGetPr(prhs[5]);
  sig = ptr[0];
  ptr = mxGetPr(prhs[6]);
  dcay1 = ptr[0];
  ptr = mxGetPr(prhs[7]);
  dcay2 = ptr[0];
  ptr = mxGetPr(prhs[8]);
  nstep = ptr[0];
  ptr = mxGetPr(prhs[9]);
  idum = (long) ptr[0];   /* seed for random generator (passed from matlab) */
  
  rlen = slen*nstep;
  /* Create memory for output variables and assign them to proper names */
  plhs[0] = mxCreateDoubleMatrix(rlen,1,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(rlen,1,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(rlen,1,mxREAL);
  Sp = mxGetPr(plhs[0]);
  V = mxGetPr(plhs[1]);
  I = mxGetPr(plhs[2]);
  
  dtsim = 1/((double) nstep);
  sig = sig*sqrt(dtsim);
  nse = gasdev(&idum);
  
  Iprev = 0;
  if (vleak < 0) 
    vprev = vleak;
  else
    vprev = 0; 
  
  for(ii=0;ii<slen;ii++)  /* Outer loop: 1 iter per time bin of input*/
    { 
      for(jj=0;jj<nstep;jj++) /* Outer loop: 1 iter per bin of output */
	{
	  kk = ii*nstep+jj; /* time index for current step */
	  II = (I[kk] + (Iprev+(Iinj[ii]-Iprev)*dtsim*(jj+1))*dtsim)*dcay2; 
	  	  
	  if (vprev > vthr) /* Spike! */
	    { 
	      II += ih[0]*dcay2;  /* update II for this bin */
	      if (kk+nh < rlen)   /* inject full h current */
		{
		  for (j=0;j<nh;j++)
		    I[kk+j] += ih[j];
		}
	      else
		{ 
		  for (j=kk;j<rlen;j++) /* inject h to only to rlen */
		    I[j] += ih[j-kk];
		}
	      /* record spike */
	      Sp[kk-1] = 1;
	      nse = gasdev(&idum);
	      V[kk] = vr*dcay1 + vleak*(1-dcay1) + II + nse*sig;
	    }
	  else
	    {
	      /* no spike */
	      nse = gasdev(&idum);
	      V[kk] = (vprev-vleak)*dcay1 + vleak + II + nse*sig;
	    }
	  vprev = V[kk];
	}
      Iprev = Iinj[ii];
    }
  if (V[rlen-1]>vthr)  /* Check for spike in last step of this time bin */
    Sp[kk] = 1;
}

