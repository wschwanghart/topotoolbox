// file: tsort_mex.c
// steepest gradient
// inp 1: SE (steepest downward neighbor grid, uint32)
//        matrix with, e.g.
//        [4 0 4 7;
//         4 4 4 7;
//         5 5 8 8]; 
//        indicates that cell 1 drains to 4, cell 2 drains 
//        to 4 and cell 9 to 8 etc... Cells with SE == 0 don't 
//        drain anywhere
// inp 2: seed pixels (linear index, uint32)
//        above example 4
// inp 3: nr of missing values (scalar, double)
//        above example 0
//
// outp 1: topological order list of linear indices (uint32)
//         above example: 12 9 6 3 11 8 5 2 10 7 1
// outp 2: linear indices of the childs of outp 1 (uint32)
//         above example: 8  8 5 5 7 7 4 4 4 4 4


#include "mex.h"
#include "matrix.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mwSize nrc, nrow, ncol; 
	mwSignedIndex IX, IXN, neigh, row, col, r, c, nrLM, pix1, pix2, pLM; //mwIndex
	mwSignedIndex coln[8] = {0, 1, 1, 1, 0, -1, -1, -1}; //int
	mwSignedIndex rown[8] = {-1, -1, 0, 1, 1, 1, 0, -1}; //int
	uint32_T *SN, *ix, *ixc, *LOCALMINIMA, *MV;

    
    if( nrhs > 3 || nlhs > 2 ) {
        mexErrMsgTxt("Need 3 inputs and no more than 2 output");
    }
    
	// input data
    SN    = mxGetData(prhs[0]);
	nrow  = mxGetM(prhs[0]);
	ncol  = mxGetN(prhs[0]);
	nrc   = nrow*ncol;

	LOCALMINIMA  = mxGetData(prhs[1]);
	nrLM  = mxGetNumberOfElements(prhs[1]);
	
	MV = mxGetData(prhs[2]);
	
	// output 
	plhs[0] = mxCreateNumericMatrix((mwSize) nrow*ncol - nrLM,(mwSize) 1, mxUINT32_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix((mwSize) nrow*ncol - nrLM,(mwSize) 1, mxUINT32_CLASS, mxREAL);
	
	ix      = (uint32_T*) mxGetData(plhs[0]);
	ixc     = (uint32_T*) mxGetData(plhs[1]);
	 
	// set initial pointer locations
	pix1 = (nrc-nrLM-MV[0]-1);
	pix2 = (nrc-nrLM-MV[0]-1);
	pLM  = 0;
	
	// IX = LOCALMINIMA(pLM)
	IX = LOCALMINIMA[pLM] - 1;
	
	while (pix1 >= 0) {
		row = IX%nrow;
		col = (IX-row)/nrow;
		
		for ( neigh = 0; neigh < 8; neigh++ ) {
			// row index of neighbor
			r = row + rown[neigh];
			// column index of neighbor
			c = col + coln[neigh];
            
			// potential neigbhor outside grid boundaries?
			// if yes, go back to beginning of the for loop
            if (r<0 || r >= nrow || c < 0 || c >=ncol) {
               continue;
            }
			
			IXN = c*nrow + r;
			
			if (SN[IXN] == (IX+1)) {
				ixc[pix1] = (uint32_T) IX+1; //int
				ix[pix1]  = (uint32_T) IXN+1; //int
				pix1 = pix1-1;
				
			}
		}
		
		// no upslope neighbor
		if (ix[pix2] == 0) {
			pLM = pLM +1;
			if (pLM > nrLM) {
			return;
			}
			IX = LOCALMINIMA[pLM] -1;
		} else {
			IX = ix[pix2] - 1;
			pix2 = pix2-1;
		}
	}
	return;


}
