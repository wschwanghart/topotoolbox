function [x,y,z] = randomsample(DEM,n)

%RANDOMSAMPLE Uniform random sampling of a GRIDobj
%
% Syntax
%
%     [x,y] = randomsample(DEM,n)
%     [x,y,z] = ...
%
% Description
%
%     randomsample creates spatially uniformly distributed random samples
%     of the surface in DEM. The function will not sample NaN-values. 
%
% Input arguments
%
%     DEM    GRIDobj
%     n      number of random samples
%
% Output arguments
%
%     x,y    nx1 coordinate vectors
%     z      nx1 vector with values obtained from DEM by linear
%            interpolation
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);     
%     D = drainagebasins(FD,S);
%     DEM.Z(D.Z==0) = nan; 
%     [x,y] = randomsample(DEM,100);
%     imageschs(DEM)
%     hold on
%     scatter(x,y,'ok','MarkerFaceColor','w','MarkerFaceAlpha',.5)
%     hold off
%
% See also: rand, GRIDobj
%        
% Author:  Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 12. July, 2018 

if nargin == 1
    n = 100;
end

I       = ~isnan(DEM.Z);
nvalid  = nnz(I);


ix    = rand(n,1)*(nvalid-1) + 1;

IX    = (uint32(1):uint32(prod(DEM.size)))';
IX    = IX(I(:));

IX = IX(round(ix));
[x,y]  = ind2coord(DEM,IX);

x = x+rand(n,1)*DEM.cellsize - DEM.cellsize/2;
y = y+rand(n,1)*DEM.cellsize - DEM.cellsize/2;

if nargout == 3
    z = interp(DEM,x,y,'linear');
end