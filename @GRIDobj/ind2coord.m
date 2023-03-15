function [x,y] = ind2coord(DEM,ix)

%IND2COORD convert linear index to x and y coordinates
%
% Syntax
%
%     [x,y] = ind2coord(DEM,ix)
%
% Description
%
%     ind2coord converts a linear index into an instance of GRIDobj to x
%     and y coordinates.
%
% Input arguments
%
%     DEM     instance of GRIDobj
%     ix      linear index
%
% Output arguments
%
%     x,y     x- and y-coordinates  
%
% See also: GRIDobj/coord2ind, GRIDobj/getcoordinates
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 30. January, 2013


[r,c] = ind2sub(DEM.size,ix(:));
xy    = double([r c ones(numel(ix),1)])*DEM.refmat;
x = xy(:,1);
y = xy(:,2);