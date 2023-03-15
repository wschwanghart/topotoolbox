function DZ = vertdistance2stream(FD,S,DEM)

%VERTDISTANCE2STREAM vertical distance to streams 
%
% Syntax
%
%     DZ = vertdistance2stream(FD,S,DEM)
%
% Description
%
%     vertdistance2stream calculates the height of each cell in a digital
%     elevation model DEM above the nearest stream cell in S along the flow
%     paths in FD (height above nearest drainage (HAND)).
%
% Input arguments
%
%     FD    instance of FLOWobj
%     DEM   digital elevation model (class: GRIDobj)
%     S     stream network (class: STREAMobj)
%
% Output arguments
%
%     DZ    vertical distance to streams (class: GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1e6,'unit','m');
%     DZ = vertdistance2stream(FD,S,DEM);
%     imageschs(DEM,DZ)
%     hold on
%     plot(S,'k','LineWidth',2)
% 
%
% See also: FLOWobj, FLOWobj/flowdistance, FLOWobj/mapfromnal, GRIDobj, 
%           STREAMobj
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 3. November, 2021


% 4/3/2016: the function now makes copies of FD.ix and FD.ixc (see 
% FLOWobj/flowacc

% 11/3/2021: performance improvement by not using GRIDobj when going
% through the loop

narginchk(3,3)

validatealignment(S,DEM);
validatealignment(FD,DEM);

Z = -inf(DEM.size,'like',DEM.Z);
Z(S.IXgrid) = DEM.Z(S.IXgrid);

ix  = FD.ix;
ixc = FD.ixc;
for r = numel(ix):-1:1
    Z(ix(r)) = max(Z(ix(r)),Z(ixc(r)));
end
DZ = DEM-Z;
DZ.name = 'Heigt above nearest drainage';
