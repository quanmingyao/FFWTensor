#include <math.h>
#include "mex.h"

/* C = spmultic_vec(I,J, A, B) computes a sparse representation(indicated by 
 subs (I,J) ) of C(:,s) where C = A(:, s) * B(s, :) */ 
#define I prhs[0]
#define J prhs[1]
#define A prhs[2]
#define B prhs[3]

#define C plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    int row, col;
    int i, j;
    double tmp;
    
    size_t length = mxGetNumberOfElements(I); 
    size_t m = mxGetM(A);
    size_t k = mxGetN(A);
    size_t n = mxGetN(B);
    
    double *Cp;
    double *Ap = mxGetPr(A);
    double *Bp = mxGetPr(B);
    
    double *Ip = mxGetPr(I);
    double *Jp = mxGetPr(J);
    
    C = mxCreateDoubleMatrix(length, k, mxREAL);
    Cp = mxGetPr(C);
    
    for (i = 0; i < length; i++) {
        row = (int)Ip[i];
        col = (int)Jp[i];
        
        for (j = 0; j < k; j++) {
            Cp[i + length * j] = Ap[(row-1) + j*m] * Bp[j + (col-1) * k];
        }
    }
}
