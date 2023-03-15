function [r,c,res] = coord2sub(DEM,x,y)

%COORD2SUB convert x and y coordinates to subscripts into a GRIDobj
%
% Syntax
%
%     [r,c] = coord2sub(DEM,x,y)
%     [r,c,res] = ...
%
% Description
%
%     coord2sub converts vectors of x and y coordinates to subscripts 
%     into an instance of GRIDobj.
%
% Input arguments
%
%     DEM     instance of GRIDobj
%     x,y     x- and y-coordinates  
%
% Output arguments
%
%     r,c     row and column subscripts
%     res     residual distance between coordinates and nearest grid cell
%             centers for coordinate pair x, y
%
% See also: GRIDobj/coord2ind, GRIDobj/sub2coord, GRIDobj/getcoordinates
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017



narginchk(3,3);
% get coordinate vectors
[X,Y] = refmat2XY(DEM.refmat,DEM.size);
X = X(:);
Y = Y(:);

% check input
np = numel(x);
if np ~= numel(y);
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

if any(I(:));
    warning('TopoToolbox:outsidegrid',...
        'There are some points outside the grid''s borders');
end

x(I)    = [];
y(I)    = [];

I = ~I;

if any(I)
    r = nan(np,1);
    r(I) = IX2(I);
    c = nan(np,1);
    c(I) = IX1(I);
    if nargout == 3
        res = nan(np,1);
        res(I) = hypot(X(IX1(I))-x,Y(IX2(I))-y);
    end
else
    r = nan(np,1);
    c = nan(np,1);
    res = nan(np,1);
end
