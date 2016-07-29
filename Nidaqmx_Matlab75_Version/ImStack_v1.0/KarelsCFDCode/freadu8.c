/* $Revision: 1.2 $ */
/* Copyright (c) 1984-98 by The MathWorks, Inc. */

/*
 * MEX-function FREADU8
 *
 * A = freadu8(fname)
 *   = freadu8(fname, offset)
 *   = freadu8(fname, offset, numBytes)
 *
 * Read bytes from the file fname into the uint8 array A.
 * Read beginning from offset (0 if not specified).
 * Read numBytes bytes (to end of file if not specified).
 *
 * Return A, a uint8 column vector containing the bytes read.
 * LENGTH(A) will give the number of bytes actually read from the file.
 *
 */

static char rcsid[] = "$Id: freadu8.c,v 1.2 1997/11/21 23:35:57 moler Exp $";

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
 *          offset   --- requested start offset for read operation
 *                       Defaults to 0 if not specified
 *          numBytes --- requested number of bytes to read
 *                       Defaults to file_length - offset if not specified
 *
 * Return:  none
 */
void ValidateInput(int nlhs, int nrhs, const mxArray *prhs[],
                   FILE **fp, int32_T *offset, int32_T *numBytes)
{
    char *fname = NULL;
    long endOfFileOffset;

    if (nrhs < 1)
    {
        mexErrMsgTxt("Too few input arguments");
    }
    if (nrhs > 3)
    {
        mexErrMsgTxt("Too many input arguments");
    }
    if (nlhs > 1)
    {
        mexErrMsgTxt("Too many output arguments");
    }

    if (!mxIsChar(prhs[0]))
    {
        mexErrMsgTxt("First input argument must be a char array");
    }
    if ((nrhs >= 2) && !mxIsDouble(prhs[1]))
    {
        mexErrMsgTxt("Second input argument must be a double scalar");
    }
    if ((nrhs >= 3) && !mxIsDouble(prhs[2]))
    {
        mexErrMsgTxt("Third input argument must be a double scalar");
    }

    if ((nrhs >= 2) && mxIsEmpty(prhs[1]))
    {
        mexErrMsgTxt("Second input argument cannot be empty");
    }

    if ((nrhs >= 3) && mxIsEmpty(prhs[2]))
    {
        mexErrMsgTxt("Third input argument cannot be empty");
    }

    fname = mxCalloc(mxGetM(prhs[0]) * mxGetN(prhs[0]) + 1, sizeof(*fname));
    mxGetString(prhs[0], fname, mxGetM(prhs[0]) * mxGetN(prhs[0]) + 1);
    
    *fp = fopen(fname, "rb");
    if (*fp == NULL)
    {
        mexErrMsgTxt(strerror(errno));
    }

    if (nrhs < 2)
    {
        *offset = 0;
    }
    else
    {
        *offset = (int32_T) mxGetScalar(prhs[1]);
    }

    if (nrhs < 3)
    {
        if (fseek(*fp, 0, SEEK_END) != 0)
        {
            fclose(*fp);
            mexErrMsgTxt(strerror(errno));
        }
        else
        {
            endOfFileOffset = ftell(*fp);
            if (endOfFileOffset == -1L)
            {
                fclose(*fp);
                mexErrMsgTxt(strerror(errno));
            }
            *numBytes = endOfFileOffset - *offset;
        }
    }
    else
    {
        *numBytes = (int32_T) mxGetScalar(prhs[2]);
    }

    mxFree((void *) fname);
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    int32_T offset;
    int32_T numBytes;
    FILE *fp;
    mxArray *A;
    mxArray *Atmp;
    mxArray *count;
    int outputSize[2];
    uint8_T *prA;
    uint8_T *prAtmp;
    size_t readCount = 0;
    int k;

    ValidateInput(nlhs, nrhs, prhs, &fp, &offset, &numBytes);

    /* Initialize output */
    outputSize[0] = numBytes;
    outputSize[1] = 1;
    A = mxCreateNumericArray(2, outputSize, mxUINT8_CLASS, mxREAL);
    prA = (uint8_T *) mxGetPr(A);
    
    /* Seek to starting offset */
    if (fseek(fp, offset, SEEK_SET) != 0)
    {
        fclose(fp);
        mexErrMsgTxt(strerror(errno));
    }

    /* Read data */
    readCount = fread((void *) prA, sizeof(uint8_T), (size_t) numBytes, fp);
    if (ferror(fp))
    {
        fclose(fp);
        mexErrMsgTxt(strerror(errno));
    }

    /* Did we read as much data as requested? */
    if (readCount < numBytes)
    {
        /* We must reallocate the output matrix with the correct size */
        outputSize[0] = readCount;
        outputSize[1] = 1;
        Atmp = mxCreateNumericArray(2, outputSize, mxUINT8_CLASS, mxREAL);
        prAtmp = (uint8_T *) mxGetPr(Atmp);
        for (k = 0; k < readCount; k++)
        {
            *prAtmp = *prA;
            prAtmp++;
            prA++;
        }
        mxDestroyArray(A);
        A = Atmp;
        Atmp = NULL;
    }

    /* Return the first output argument and clean up */
    plhs[0] = A;

    fclose(fp);
}
