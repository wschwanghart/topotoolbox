function [x,y] = getcoordinates(DEM,type)

%GETCOORDINATES get coordinate vectors of an instance of GRIDobj
%
% Syntax
%
%     [x,y] = getcoordinates(DEM)
%     [x,y] = getcoordinates(DEM,type)
%
% Input arguments
%
%     DEM    grid (class: GRIDobj)
%
% Output arguments
%
%     x      coordinate vector in x direction (row vector)
%     y      coordinate vector in y direction (column vector)
%     type   'vector' (default). Alternatively, you can return coordinate
%            matrices ('matrix') or GRIDobjs ('GRIDobj')
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     [x,y] = getcoordinates(DEM);
%     surf(x,y,double(DEM.Z))
%     axis image; shading interp; camlight
%     
%
%
% See also: GRIDobj2mat
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. September, 2018

[x,y] = refmat2XY(DEM.refmat,DEM.size);

if nargin == 1
    return;
else
    type = validatestring(type,{'vector','matrix','GRIDobj'});
end

switch type
    case 'vector'
        return
    case 'matrix'
        [x,y] = meshgrid(x,y);
    case 'GRIDobj'
        [x,y] = meshgrid(x,y);
        X = GRIDobj(DEM);
        X.Z = x;
        Y = GRIDobj(DEM);
        Y.Z = y;
        x = X;
        y = Y;
end