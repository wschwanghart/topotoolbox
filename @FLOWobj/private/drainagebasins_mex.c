// file: drainagebasins_mex.c
// drainage basin delineation using the FlowDirObj class
#include "mex.h"
#include "matrix.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mwSize numel; // use mwIndex
	mwIndex index, ix, ixc;
    uint32_T *ixp;
	uint32_T *ixcp;
	uint32_T DBcounter;
	uint32_T *D;
	double *siz;
    
    if( nrhs > 3 || nlhs > 1 ) {
        mexErrMsgTxt("Need 3 inputs and no more than 1 output");
    }
    
    ixp    = mxGetData(prhs[0]);
	ixcp   = mxGetData(prhs[1]);
	siz    = mxGetPr(prhs[2]);
	numel  = mxGetNumberOfElements(prhs[0]);
	
	/* Create an m-by-n mxArray and either copy initial values into it or create a new one with zeros*/
    // double precision
	// plhs[0] = mxCreateDoubleMatrix((mwSize) siz[0],(mwSize) siz[1], mxREAL);
    /* pointer into drainage basins*/
	// D       = mxGetPr(plhs[0]);
    
    // uint32
    plhs[0] = mxCreateNumericMatrix((mwSize) siz[0],(mwSize) siz[1], mxUINT32_CLASS, mxREAL);
    D       = mxGetData(plhs[0]);
	
	/* set pointer to last element in index arrays */
	ixp  = ixp + numel;
	ixcp = ixcp + numel;
	DBcounter = 0;
	
    for( index= 0; index<=numel; index++ ) {
		
		ix  = (mwIndex) *ixp--  -1;
		ixc = (mwIndex) *ixcp-- -1;
		
		if (D[ixc] == 0) {
			DBcounter = DBcounter + 1;
			D[ixc] = DBcounter;
		}
		
		D[ix] = D[ixc];
    }
    
    return;
}