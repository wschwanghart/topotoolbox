// file: gradient8_mex.c
// steepest gradient
// inp 1: DEM
// inp 2: GWDT (output of grayweighted distance transform)
//        nonflat areas are inf
//        sills are zero
// inp 3: cellsize (scalar)
//
// outp 1: index of steepest neighbor

#include "mex.h"
#include "matrix.h"
#include "math.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mwSize nrow, ncol; // use mwIndex
	mwSignedIndex ix, neigh, row, col, r, c, ixn; //mwIndex
	mwSignedIndex coln[8] = {0, 1, 1, 1, 0, -1, -1, -1}; //int
	mwSignedIndex rown[8] = {-1, -1, 0, 1, 1, 1, 0, -1}; //int
	double z, zn, g, G, ggwdt, Ggwdt;
	double *cs;
	double *DEM, *GWDT;
	uint32_T *SN; //int
	double dist[8] = {1, sqrt(2), 1, sqrt(2), 1, sqrt(2), 1, sqrt(2) };

    
    if( nrhs > 3 || nlhs > 1 ) {
        mexErrMsgTxt("Need 3 inputs and no more than 1 output");
    }
    
    DEM    = mxGetPr(prhs[0]);
    GWDT   = mxGetPr(prhs[1]);
	nrow   = mxGetM(prhs[0]);
	ncol   = mxGetN(prhs[0]);
	cs     = mxGetPr(prhs[2]);
	
	/* Create an m-by-n mxArray and either copy initial values into it or create a new one with ones*/
	plhs[0] = mxCreateNumericMatrix((mwSize) nrow,(mwSize) ncol, mxUINT32_CLASS, mxREAL);
	/* pointer into drainage basins*/
	SN      = mxGetData(plhs[0]);
	
	// set cell distances
	for( neigh = 0; neigh < 8; neigh++) {
		dist[neigh] = cs[0] * dist[neigh];
	}

    for( row= 0; row<nrow; row++ ) {
	
		for( col= 0; col<ncol; col++) {
			
			ix = (col*nrow) + row;
			z  = DEM[ix];
			G = 0;
            Ggwdt = GWDT[ix];
            
			// nan?
			if (mxIsNaN(z)) {
				SN[ix] = 0;
                continue;
			}
		
			for ( neigh = 0; neigh < 8; neigh++ ) {
				// row index of neighbor
				r = row + rown[neigh];
				// column index of neighbor
				c = col + coln[neigh];
                
                if (r<0 || r >= nrow || c < 0 || c >=ncol) {
                   continue;
                }
				
				ixn = (c*nrow) + r;
				// calculate slope
				g = (z - DEM[ixn])/ dist[neigh];
				// find maximum slope
				if ( g > G ) {
				   SN[ix] = ixn + 1;
				   G = g;
				} else if ( g == 0 ) {
                   // if slope is 0
                   if ( mxIsInf(GWDT[ix]) || GWDT[ix] == 0) {
                       // we are not in a flat section
                       continue;
                   } else {
                       // we are in a flat section
                       if (GWDT[ixn] == 0) {
                           // the neighbor is the sill
                           SN[ix] = ixn + 1;
                           break;
                       } else {
                           ggwdt = GWDT[ixn];
                           if ( ggwdt < Ggwdt) {
                               SN[ix] = ixn + 1;
                               Ggwdt = ggwdt;
                           }
                       }
                   }
                }
			}
		}        
    }
    return;
}
