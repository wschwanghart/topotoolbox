function DEM = erode(DEM,SE)

%ERODE morphological erosion
%
% Syntax
%
%     DEMd = erode(DEM,SE)
%
% Description
%
%     erode is a simple wrapper around imerode that handles nans.
%
% See also: IMDILATE
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 6. February, 2013


if isa(DEM.Z,'double') || isa(DEM.Z,'single')
    I = isnan(DEM.Z);
    DEM.Z(I) = inf;
    checknans = true;
else
    checknans = false;
end

DEM.Z = imerode(DEM.Z,SE);
if checknans
    DEM.Z(I) = nan;
end
