function G = FLOWobj2GRIDobj(F)

%FLOWOBJ2GRIDOBJ create ESRI ArcGIS flow direction grid from FLOWobj
%
% Syntax
%
%     AF = FLOWobj2GRIDobj(FD)
%
% Description
%
%     FLOWobj2GRIDobj converts an instance of FLOWobj to a grid format used
%     by ESRI ArcGIS to store flow directions. The function should be used
%     together with GRIDobj2geotiff if you wish to work in ArcGIS with flow
%     directions derived with TopoToolbox.
%     
% Input arguments
%
%     FD    instance of FLOWobj
%
% Output arguments
%
%     AF    instance of GRIDobj that stores the flow directions in the
%           format used by ESRI ArcGIS.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     ArcFlow = FLOWobj2GRIDobj(FD);
%     imagesc(ArcFlow)
%
% 
% See also: GRIDobj, GRIDobj2geotiff
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. October, 2013

switch lower(F.type)
    case {'multi' 'dinf'}
        error('not possible for multiple flow directions')
end

G = copy2GRIDobj(F);

% ESRI's flow direction raster is defined here:
% http://help.arcgis.com/de/arcgisdesktop/10.0/help/index.html#//009z00000052000000
%
% The directions are 
%
%  ESRI             FLOWobj
%
%  32 64 128        -1-nrrows -1 -1+nrrows
%  16     1     =   -nrrows    0   nrrows
%   8  4  2         1-nrrows   1  1+nrrows

nrrows = G.size(1);
ix     = int32(F.ix);
ixc    = int32(F.ixc);
direction = int32([nrrows 1+nrrows 1 1-nrrows -nrrows -1-nrrows -1 -1+nrrows]');
[~,fdESRI] = ismember(ixc-ix,direction);
fdESRI = uint8(2.^(fdESRI-1));

G.Z = zeros(G.size,'uint8');
G.name = 'flow direction grid';
G.Z(F.ix) = fdESRI;

% The definition of flow direction for cells with internal drainage is
% somewhat strange
% ESRI help on flow direction:
% If a cell is lower than its eight neighbors, that cell is given the value
% of its lowest neighbor, and flow is defined toward this cell. If multiple
% neighbors have the lowest value, the cell is still given this value, but
% flow is defined with one of the two methods explained below. This is used
% to filter out one-cell sinks, which are considered noise.

% Here we mimic this behavior by reversing the flow direction towards one
% of the neighbors
I   = G.Z(ixc)==0;
ixc = ixc(I);
ix  = ix(I);
[~,fdESRI] = ismember(ix-ixc,direction);
fdESRI = uint8(2.^(fdESRI-1));
[ixc,ia] = unique(ixc);
fdESRI = fdESRI(ia);
G.Z(ixc) = fdESRI;






