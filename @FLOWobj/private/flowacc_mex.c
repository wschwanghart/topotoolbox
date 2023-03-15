// file: flowacc_mex.c
// flow accumulation using the FlowDirObj class
#include "mex.h"
#include "matrix.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mwSize numel, numelW0, ncol, nrow; // use mwIndex
	mwIndex index, ix, ixc;
    uint32_T *ixp;
	uint32_T *ixcp;
	double *Ap;
	double *W0p;
    
    if( nrhs > 3 || nlhs > 1 ) {
        mexErrMsgTxt("Need 3 inputs and no more than 1 output");
    }
    
    ixp    = (uint32_T*) mxGetData(prhs[0]); //int*
	ixcp   = (uint32_T*) mxGetData(prhs[1]); //int*
	W0p    = mxGetPr(prhs[2]);
	numel    = mxGetNumberOfElements(prhs[0]);
	numelW0  = mxGetNumberOfElements(prhs[2]);
	
	/* Create an m-by-n mxArray and either copy initial values into it or create a new one with ones*/
	if(numelW0==2) {

		plhs[0] = mxCreateDoubleMatrix((mwSize) W0p[0],(mwSize) W0p[1], mxREAL);
		Ap = mxGetPr(plhs[0]);
		for(index = 0; index<(W0p[0]*W0p[1]); index++) {
			Ap[index] = 1;
		}
	} else {
		ncol   = mxGetN(prhs[2]);
		nrow   = mxGetM(prhs[2]);
		plhs[0] = mxCreateDoubleMatrix(nrow,ncol, mxREAL);
		Ap = mxGetPr(plhs[0]);
		for(index = 0; index<numelW0; index++) {
			Ap[index] = W0p[index];
		}
	}
	

    for( index=0; index<numel; index++ ) {
		ix  = (mwIndex) *ixp++ -1;
		ixc = (mwIndex) *ixcp++ -1 ;
		Ap[ixc] = Ap[ix] + Ap[ixc];
    }
    return;
}