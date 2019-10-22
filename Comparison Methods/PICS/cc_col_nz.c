/*
 * cc_col_nz.c  Faster MEX implementation of cc_col_nz.m
 * Run "mex cc_col_nz.c" on mcfarland (NEED Correct LD_LIBRARY_PATH !)
 *
 */

/* $Id: cc_col_nz.c,v 1.3 2004/04/28 19:23:53 deepay Exp $ */

/*
 * $Log: cc_col_nz.c,v $
 * Revision 1.3  2004/04/28 19:23:53  deepay
 * stuff for mex
 *
 * Revision 1.2  2004/02/20 19:51:21  spapadim
 * Add bounds check for cluster labels (fixes coredump upon invalid calls)
 *
 * Revision 1.1  2004/02/20 10:22:12  spapadim
 * Initial checkin
 *
 */

#include "mex.h"

#define C_IN   prhs[0]
#define k_IN   prhs[1]
#define Qx_IN  prhs[2]
#define nz_OUT plhs[0]

#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define MIN(a,b)  (((a) < (b)) ? (a) : (b))

void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *Qx, *nz;
    int *C_irow, *C_jcol;
    int k, C_nnz;
    unsigned int i;

    if (nrhs != 3)
        mexErrMsgTxt("Wrong number of input arguments");
    if (nlhs > 1)
        mexErrMsgTxt("Too many output arguments");
    if (!mxIsSparse(C_IN))
        mexErrMsgTxt("C must be sparse");
    if (mxGetN(C_IN) != 1)
        mexErrMsgTxt("C must be a sparse column vector");
    
    k = (int)mxGetScalar(k_IN);
    Qx = mxGetPr(Qx_IN);
    C_irow = mxGetIr(C_IN);
    C_jcol = mxGetJc(C_IN);
    /* C_nzmax = mxGetNzmax(C_IN); */
    C_nnz = C_jcol[1];

    /* mexPrintf("nzmax %d, nnz %d\n", C_nzmax, C_nnz); */

    nz_OUT = mxCreateDoubleMatrix(k, 1, mxREAL);
    nz = mxGetPr(nz_OUT);

    for (i = 0;  i < C_nnz; i++) {
        int row_idx = (int)C_irow[i];
        int cluster_label = (int)Qx[row_idx];
        if ((cluster_label < 1) || (cluster_label > k))
            mexErrMsgTxt("Invalid cluster label value");
        /* mexPrintf("Row idx: %d\n", row_idx); */
        nz[cluster_label-1] = nz[cluster_label-1] + 1.0;
    }
}
