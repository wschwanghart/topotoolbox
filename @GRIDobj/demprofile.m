function [dn,z,x,y] = demprofile(DEM,n,x,y)

%DEMPROFILE get profile along path
%
% Syntax
%
%     [d,z] = demprofile(DEM)
%     [d,z] = demprofile(DEM,n)
%     [d,z,x,y] = demprofile(DEM,n,x,y)
%     p     = demprofile(DEM,...)
%     p     = demprofile(DEM,p)
%
% Description
%
%     demprofile enables to interactively or programmetically derive
%     elevation profiles from DEMs.
%
% Input arguments
%
%     DEM    digital elevation model (class: GRIDobj)
%     n      number of points along profile
%     x,y    coordinate vectors of profile vertices
%     p      structure array (struct) returned by demprofile
%
% Output arguments
%
%     d      distance along profile
%     z      elevation values interpolated (linear) onto profile
%     x,y    coordinate vectors of profile 
%
%     p      If only one output argument is used than demprofile returns
%            a scalar structure array with fieldnames that refer to above 
%            listed output arguments.
% 
% See also: GRIDobj/interp, GRIDobj/measure
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 2. September, 2020


% interactive
if nargin <= 2
    
    if nargin == 2 && isstruct(n)
        x = n.x;
        y = n.y;
        dn = getdistance(x,y);
    else
    
    h = imagesc(DEM);
    ax = get(h,'Parent');

    [x, y] = getline(ax);
    do = getdistance(x,y);
    
    if nargin == 1
        n = ceil(do(end)/DEM.cellsize)*2;
    end
    
    dn = linspace(0,do(end),n);
    dn = dn(:);
    xy = interp1(do,[x y],dn,'linear');
    
    x  = xy(:,1);
    y  = xy(:,2);
    end
    
else
    if n ~= numel(x);
        x = x(:);
        y = y(:);
        do = getdistance(x,y);
        dn = linspace(0,do(end),n);
        dn = dn(:);
        xy = interp1(do,[x y],dn,'linear');
    
        x  = xy(:,1);
        y  = xy(:,2);
    else
    
        dn = getdistance(x,y);
    end
end

z = double(interp(DEM,x,y));

if nargout == 1
    
    out.d = dn;
    out.z = z;
    out.x = x;
    out.y = y;
    
    dn = out;
end

end
    


function cumdxy = getdistance(x,y)
    x = x(:);
    y = y(:);

    dx = diff(x);
    dy = diff(y);
    dxy = hypot(dx, dy); % square root of sum of squares
    cumdxy = [0; cumsum(dxy,1)];
end % 