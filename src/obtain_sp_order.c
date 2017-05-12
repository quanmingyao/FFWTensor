#include <math.h>
#include "mex.h"

#define A prhs[0]
#define IDX plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mwSize nzmax = mxGetNzmax(A);
    double *A_Pr;
    double *IDX_Pr;
    
    IDX = mxCreateDoubleMatrix(nzmax, 1, mxREAL);
    IDX_Pr = mxGetPr(IDX);

    A_Pr = mxGetPr(A);  

    memcpy(IDX_Pr, A_Pr, (size_t)(sizeof(mwIndex)*nzmax));
}
