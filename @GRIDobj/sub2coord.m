function [x,y] = sub2coord(DEM,r,c)

%SUB2COORD convert subscripts to x and y coordinates
%
% Syntax
%
%     [x,y] = sub2coord(DEM,r,c)
%
% Description
%
%     sub2coord converts a subscripts into an instance of GRIDobj to x
%     and y coordinates.
%
% Input arguments
%
%     DEM     instance of GRIDobj
%     r,c     row and column indices
%
% Output arguments
%
%     x,y     x- and y-coordinates  
%
% See also: GRIDobj/coord2sub, GRIDobj/getcoordinates
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 30. January, 2013


xy    = double([r c ones(numel(r),1)])*DEM.refmat;
x = xy(:,1);
y = xy(:,2);