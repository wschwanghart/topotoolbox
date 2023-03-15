function zi = interp(DEM,xi,yi,method)

%INTERP interpolate to query locations
%
% Syntax
%
%     zi = interp(DEM,xi,yi)
%     zi = interp(DEM,xi,yi,method)
%
% Description
%
%     INTERP uses the griddedInterpolant class to interpolate values in the
%     instance of GRIDobj (DEM) to query locations at xi and yi. 
%     If DEM.Z is an integer class, interp will convert it to single
%     precision and use linear interpolation as default. If DEM.Z is
%     logical, nearest neigbhor will be used by default.
%
% Input arguments
%
%     DEM     instance of GRIDobj
%     xi,yi   x- and y-coordinates of query locations 
%     method  interpolation method (default = 'linear'). See the
%             documentation of the griddedInterpolant class for further
%             methods
% 
% Output arguments
%
%     zi      interpolated values at query locations
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     [x,y] = getoutline(DEM);
%     xy = rand(20,2);
%     xy(:,1) = xy(:,1)*(max(x)-min(x)) + min(x);
%     xy(:,2) = xy(:,2)*(max(y)-min(y)) + min(y);
%     z = interp(DEM,xy(:,1),xy(:,2));
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017


narginchk(3,4)

if nargin == 3
    method = 'linear';
else
    method = validatestring(method,...
        {'linear','nearest','spline','pchip','cubic'},'interp','method',4);
end

% created griddedInterpolant class
[x,y] = getcoordinates(DEM);
% flip y to get monotonic increasing grid vectors
y     = y(end:-1:1);

if isinteger(DEM.Z)
    cla   = class(DEM.Z);
    DEM.Z = single(DEM.Z);
    convoutput = true;
elseif islogical(DEM.Z)
    cla   = class(DEM.Z);
    DEM.Z = single(DEM.Z);
    if nargin == 3;
        method = 'nearest';
    end
    convoutput = true;
else
    convoutput = false;
end


F     = griddedInterpolant({x,y},flipud(DEM.Z)',method); 

% interpolate
zi     = F(xi,yi);

if convoutput
    zi = cast(zi,cla);
end


