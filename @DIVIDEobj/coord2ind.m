function [IX,res] = coord2ind(D,x,y)

%COORD2IND convert x and y coordinates to linear index
%
% Syntax
%
%     IX = coord2ind(D,x,y)
%     [IX,res] = ...
%
% Description
%
%     coord2ind converts vectors of x and y coordinates to a linear index 
%     into an instance of DIVIDEobj. Note that the virtual grid associated
%     with a DIVIDEobj is shifted by half a cell size in x and y and has an
%     extra row and collumn. The linear indices point to the coreners of
%     the pixels in the GRIDobj the DIVIDEobj is based on.
%
% Input arguments
%
%     D       instance of DIVIDEobj
%     x,y     x- and y-coordinates  
%
% Output arguments
%
%     ix      linear index
%     res     residual distance between coordinates and nearest grid cell
%             centers for coordinate pair x, y
%
% See also: GRIDobj/coord2ind, GRIDobj/ind2coord, GRIDobj/getcoordinates
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
%         code simply copied and help adjusted by Dirk Scherler
% Date: 17. August, 2017 and February 2019



narginchk(3,3);
% get coordinate vectors
[X,Y] = refmat2XY(D.refmat,D.size);
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

I = IX1>D.size(2) | IX1<1 | IX2>D.size(1) | IX2<1;

if any(I(:))
    warning('TopoToolbox:outsidegrid',...
        'There are some points outside the grid''s borders');
end

% Check for nans because they cause errors since 2024a
I = I | isnan(x) | isnan(y);

x(I)    = [];
y(I)    = [];

I = ~I;

if any(I)
    IX = nan(np,1);
    IX(I) = sub2ind(D.size,IX2(I),IX1(I));
    if nargout == 2
        res = nan(np,1);
        res(I) = hypot(X(IX1(I))-x,Y(IX2(I))-y);
    end
else
    IX = nan(np,1);
    res = nan(np,1);
end




