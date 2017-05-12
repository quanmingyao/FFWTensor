#include <math.h>
#include "mex.h"

/* C = spmultic(I,J, A, B) computes a sparse representation(indicated by 
 subs (I,J) ) of C where C = A*B */ 
#define I prhs[0]
#define J prhs[1]
#define A prhs[2]
#define B prhs[3]

#define C plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    //Obtain the nonzero entry numbers
    
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
    
    C = mxCreateDoubleMatrix(length, 1, mxREAL);
    Cp = mxGetPr(C);
    
    //printf("%d", length);
    
    
    for (i = 0; i < length; i++) {
        row = (int)Ip[i];
        col = (int)Jp[i];
        
        tmp;
        tmp = 0;
        //printf("the %d th element\n", i);
        for (j = 0; j < k; j++){
            //printf("%f multiplies %f\n", Ap[(row-1) + j*m], Bp[j + (col-1) * k]);
            tmp = tmp + Ap[(row-1) + j*m] * Bp[j + (col-1) * k];
        }
        
        Cp[i] = tmp;
    }
     
}




