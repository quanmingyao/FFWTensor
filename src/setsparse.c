#include <math.h>
#include "mex.h"

#define A prhs[0]
#define VAL prhs[1]
#define I prhs[2]
#define J prhs[3]

#define B plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    size_t m = mxGetM(A);
    size_t n = mxGetN(A);
    mwSize nzmax = mxGetNzmax(A);
    mwIndex *B_Ir, *B_Jc;
    double *B_Pr;
    mwIndex *A_Ir, *A_Jc;
    double *A_Pr;
    double *VALp;
    
    B = mxCreateSparse((mwSize)m, (mwSize)n, nzmax, mxREAL);

    B_Ir = mxGetIr(B);
    B_Jc = mxGetJc(B);
    B_Pr = mxGetPr(B);
    
    A_Ir = mxGetIr(A);
    A_Jc = mxGetJc(A);
    A_Pr = mxGetPr(A);  
    
    VALp = mxGetPr(VAL);
    memcpy(B_Ir, A_Ir, (size_t)(sizeof(mwIndex)*nzmax) );
    memcpy(B_Jc, A_Jc, (size_t)(sizeof(mwIndex)*(n+1)) );
    memcpy(B_Pr, VALp, (size_t)(sizeof(double)*nzmax) ); 
}
