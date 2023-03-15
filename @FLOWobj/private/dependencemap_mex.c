// file: dependencemap_mex.c
// dependendence map using the FlowDirObj class
//
// input arguments
// ix, ixc, siz, ixs
#include "mex.h"
#include "matrix.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mwSize numel, nixs; // use mwIndex
	mwIndex index, ix, ixc, ixs;
    uint32_T *ixp;
	uint32_T *ixcp;
	double *ixsp;
	mxLogical *D;
	double *siz;
    
    if( nrhs > 4 || nlhs > 1 ) {
        mexErrMsgTxt("Need 4 inputs and no more than 1 output");
    }
    
    ixp    = mxGetData(prhs[0]);
	ixcp   = mxGetData(prhs[1]);
	siz    = mxGetPr(prhs[2]);
	numel  = mxGetNumberOfElements(prhs[0]);
	ixsp   = mxGetPr(prhs[3]);
	nixs   = mxGetNumberOfElements(prhs[3]);
	
	/* Create an m-by-n mxArray and either copy initial values into it or create a new one with ones*/
	plhs[0] = mxCreateLogicalMatrix((mwSize) siz[0],(mwSize) siz[1]);
	/* pointer into drainage basins*/
	D       = mxGetLogicals(plhs[0]);
	for( index = 0; index < nixs; index++) {
        D[(mwIndex) *ixsp -1] = true;
        ixsp++;
		}
	
	/* set pointer to last element in index arrays */
	ixp  = ixp + numel;
	ixcp = ixcp + numel;
	
    for( index= 0; index<=numel; index++ ) {
		
		ix  = (mwIndex) *ixp--  -1;
		ixc = (mwIndex) *ixcp-- -1;
		
		if (D[ixc]) {
			D[ix] = D[ixc];
		}
    }
    return;
}