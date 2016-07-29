/* $Revision: 1.2 $ */
/* Copyright (c) 1984-98 by The MathWorks, Inc. */

/*
 * MEX-function FWRITEU8
 *
 * count = fwriteu8(fname, A)
 *
 * Append the elements of the uint8 array A to the file fname.
 * Return the count of bytes actually written.
 *
 */

static char rcsid[] = "$Id: fwriteu8.c,v 1.2 1997/11/21 23:35:59 moler Exp $";

#include <errno.h>
#include <string.h>
#include "mex.h"

/*
 * ValidateInput --- input argument parsing and error checking
 *
 * Inputs:  nlhs --- number of left-side arguments
 *          nrhs --- number of right-side arguments
 *          prhs --- array of pointers to right-side arguments
 *
 * Outputs: fp       --- file pointer
 *          A        --- data array
 *
 * Return:  none
 */
void ValidateInput(int nlhs, int nrhs, const mxArray *prhs[],
                   FILE **fp, const mxArray **A)
{
    char *fname = NULL;
    long length;

    if (nrhs < 2)
    {
        mexErrMsgTxt("Too few input arguments");
    }
    if (nrhs > 2)
    {
        mexErrMsgTxt("Too many input arguments");
    }
    if (nlhs > 1)
    {
        mexErrMsgTxt("Too many output arguments");
    }

    if (!mxIsChar(prhs[0]))
    {
        mexErrMsgTxt("FILENAME must be a string");
    }

    if (!mxIsUint8(prhs[1]))
    {
        mexErrMsgTxt("Class of data array must be uint8");
    }

    length = mxGetM(prhs[0]) * mxGetN(prhs[0]) + 1;
    fname = mxCalloc(length, sizeof(*fname));
    mxGetString(prhs[0], fname, length);
    
    *fp = fopen(fname, "ab");
    if (*fp == NULL)
    {
        mexErrMsgTxt(strerror(errno));
    }

    *A = prhs[1];

    mxFree((void *) fname);
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    size_t numBytes;
    size_t lengthA;
    FILE *fp;
    const mxArray *A;
    uint8_T *prA;

    ValidateInput(nlhs, nrhs, prhs, &fp, &A);

    lengthA = (size_t) (mxGetM(A) * mxGetN(A));
    prA = (uint8_T *) mxGetPr(A);

    numBytes = fwrite((void *) prA, sizeof(uint8_T), lengthA, fp);
    if (numBytes < lengthA)
    {
        mexWarnMsgTxt(strerror(errno));
    }

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *((double *) mxGetPr(plhs[0])) = (double) numBytes;

    fclose(fp);
}
