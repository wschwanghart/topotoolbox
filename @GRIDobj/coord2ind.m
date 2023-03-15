function [IX,res,valid] = coord2ind(DEM,x,y)

%COORD2IND convert x and y coordinates to linear index
%
% Syntax
%
%     IX = coord2ind(DEM,x,y)
%     [IX,res] = ...
%     [IX,res,valid] = ...
%
% Description
%
%     coord2ind converts vectors of x and y coordinates to a linear index 
%     into an instance of GRIDobj.
%
% Input arguments
%
%     DEM     instance of GRIDobj
%     x,y     x- and y-coordinates  
%
% Output arguments
%
%     ix      linear index
%     res     residual distance between coordinates and nearest grid cell
%             centers for coordinate pair x, y
%     valid   logical vector same size as x (or y) with true where x and y
%             pairs are inside the grid, and otherwise false. 
%
% See also: GRIDobj/ind2coord, GRIDobj/getcoordinates
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 20. February, 2020



narginchk(3,3);
% get coordinate vectors
[X,Y] = refmat2XY(DEM.refmat,DEM.size);
X = X(:);
Y = Y(:);

% check input
np = numel(x);
if np ~= numel(y)
    error('TopoToolbox:wronginput',...
        'x and y must have the same number of elements.');
end

% force column vectors
x = x(:);
y = y(:);

dx  = X(2)-X(1);
dy  = Y(2)-Y(1);

IX1 = (x-X(1))./dx + 1;
IX2 = (y-Y(1))./dy + 1;

IX1 = round(IX1);
IX2 = round(IX2);

I = IX1>DEM.size(2) | IX1<1 | IX2>DEM.size(1) | IX2<1;

if any(I(:))
    warning('TopoToolbox:outsidegrid',...
        'There are some points outside the grid''s borders');
end

x(I)    = [];
y(I)    = [];

I = ~I;

if nargout == 3
    valid = I;
end


if any(I)
    IX = nan(np,1);
    IX(I) = sub2ind(DEM.size,IX2(I),IX1(I));
    if nargout >= 2
        res = nan(np,1);
        res(I) = hypot(X(IX1(I))-x,Y(IX2(I))-y);
    end
else
    IX = nan(np,1);
    res = nan(np,1);
end




