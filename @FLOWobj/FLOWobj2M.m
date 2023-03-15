function M = FLOWobj2M(FD)

%FLOWOBJ2M convert instance of FLOWobj to flow direction matrix 
%
% Syntax
%
%     M = FLOWobj2M(FD);
%
% Description
%
%     FLOWobj2M converts an instance of FLOWobj to the flow direction
%     matrix M as used in TopoToolbox 1.
%
% Input arguments 
%
%     FD     FLOWobj
%
% Output arguments
%
%     M      sparse transfer matrix
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     M = FLOWobj2M(FD);
%     [x,y] = getcoordinates(DEM);
%     [x,y] = meshgrid(x,y);
%     gplot(M,[x(:) y(:)]) 
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. August, 2017


nrc = prod(FD.size);
switch FD.type
    case 'multi'
        M = sparse(double(FD.ix),double(FD.ixc),FD.fraction,nrc,nrc);
    case 'single'
        M = sparse(double(FD.ix),double(FD.ixc),1,nrc,nrc);
end
end