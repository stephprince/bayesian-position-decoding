

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>



void Usage(void) ;


/******************************************************************************
  INTERFACE FUNCTION
  */

void mexFunction(
    int           nlhs,           /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int           nrhs,           /* number of inputs */
    const mxArray *prhs[])         /* array of pointers to input arguments */
{

      double 			  *dblptr;
      int               i,j,k,z;
      double      *values;
      double      *lookupvalues;
      double			  *output;
      int				nvalues  , nlookupvalues, nrows, ncols, nlookuprows, nlookupcols, issame, abort;

      /* Check numbers of arguments */
      if (!(nrhs == 2) || !((nlhs == 1)||(nlhs == 0))) {
         Usage();
      }
   

      values = mxGetData(prhs[0]);      
      nrows = mxGetM(prhs[0]);
      ncols = mxGetN(prhs[0]);
      lookupvalues = mxGetData(prhs[1]);      
      nlookuprows = mxGetM(prhs[1]);
      nlookupcols = mxGetN(prhs[1]);
      
      if (!(ncols == nlookupcols)){
         abort = 1;
      }
      else {
         abort = 0;
      }  
    
 
      output = (double *) mxCalloc(nrows, sizeof(double));

      for (i = 0; i < nrows; i++) {
         if (abort == 0) {
            for (j = 0; j < nlookuprows; j++) {
               issame = 1;
               for (k = 0; k < ncols; k++) { 
                  if (values[i+(nrows*k)] != lookupvalues[j+(nlookuprows*k)]) {
                     issame = 0;
                  }
               }
               if (issame) {
                  output[i] = j+1;
                  break;
               }   
            }
            if (!issame) {
               output[i] = 0;
            }
         }
         else {
            output[i] = 0;
         }
      }

      plhs[0] = mxCreateDoubleMatrix(nrows,1, mxREAL);

      dblptr = mxGetPr(plhs[0]);

      for (i = 0; i < nrows; i++) {
               *dblptr = output[i];
               dblptr++;
      }


      mxFree(output);

      return;
}

void Usage(void)
{
      mexErrMsgTxt("Usage: one putput and two inputs\n");
}


