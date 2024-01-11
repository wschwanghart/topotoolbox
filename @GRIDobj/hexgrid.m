function [X,Y] = hexgrid(DEM,d)

%HEXGRID creates an array of haxagonal points
%
% Syntax
%
%     [X,Y] = hexgrid(DEM,d)
%
% Description
%
%     hexgrid creates an array of haxagonal points that are evenly
%     distributed across the entire DEM. The X and Y arrays are the center
%     coordinates of the grid cells. Currently only limited options.
%
% Input arguments
%
%     DEM   GRIDobj
%     d     factor (scalar) by which the grid cellsize is multiplied in the
%           x-direction to yield the desired point spacing
%
% Output arguments
%
%     X,Y   x,y coordinate arrays
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     [X,Y] = hexgrid(DEM,100);
%     figure
%     imagesc(DEM)
%     hold on
%     scatter(X(:),Y(:),50,'ko','filled')
%     [VX,VY] = voronoi(X,Y);
%     plot(VX,VY,'w-','LineWidth',1)
%     axis(getextent(DEM))
%
% See also: meshgrid
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: 18. August 2023

[x,y] = getcoordinates(DEM);
cs = DEM.cellsize;

dx = cs*d;
dy = dx*sqrt(3)/2;

minx = min(x);
maxx = max(x);
miny = min(y);
maxy = max(y);

x = minx:dx:maxx;
y = miny:dy:maxy;

[X,Y] = meshgrid(x,y);

[nr,nc] = size(X);
sx = zeros(nr,1);
sx(1:2:end) = sx(1:2:end)+1;
dX = repmat(sx.*dx/2,[1,nc]);
X = X + dX;


end