function DEM = imposemin(FD,DEM,sl)

%IMPOSEMIN minima imposition (carving) along drainage network
%
% Syntax
%
%     DEMc = imposemin(FD,DEM)
%     DEMc = imposemin(FD,DEM,sl)
%
% Description
%
%     imposemin carves a digital elevation model by downstream minima
%     imposition such that
%
%     DEM(ixc) = min(DEM(ix),DEM(ixc))
%
%     where ix is a linear index in DEM and ixc is the linear index of
%     the downward neighbor of ix. Depending on the breach length, long
%     flat channel section may be generated which can be avoided by
%     imposing a slight, minimum downward gradient sl such that
%
%     DEM(ixc) = min(DEM(ix)-dx*sl,DEM(ixc)) 
%
%     Take care to not choose a high gradient which will cause that
%     channels dip below the true land surface. 
%
% Input arguments
%
%     FD        flow direction object (FLOWObj)
%     DEM       digital elevation model (GRIDobj)
%     sl        minimum gradient [m/m] in downward direction (e.g. 0.001)
%
% Output arguments
%
%     DEMc      carved DEM
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     DEMc = imposemin(FD,DEM);
%     imageschs(DEM,DEM-DEMc,'colormap',flipud(parula)) 
%
% See also: STREAMobj/imposemin
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. August, 2017


% 4/3/2016: the function now makes copies of FD.ix and FD.ixc (see 
% FLOWobj/flowacc


narginchk(2,3)
validatealignment(FD,DEM);
dem = DEM.Z;

ix = FD.ix;
ixc = FD.ixc;

if nargin == 2;    
    for r = 1:numel(ix);
        dem(ixc(r)) = min(dem(ix(r)),dem(ixc(r)));
    end
    
elseif nargin == 3;
    d = getdistance(FD.ix,FD.ixc,FD.size,FD.cellsize);
    d = d*sl;
    for r = 1:numel(ix);
        dem(ixc(r)) = min(dem(ix(r))-d(r),dem(ixc(r)));
    end
end

DEM.Z = dem;
DEM.name = [DEM.name ' - imposemin'];


