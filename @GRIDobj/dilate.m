function DEM = dilate(DEM,SE)

%DILATE morphological dilation
%
% Syntax
%
%     DEMd = dilate(DEM,SE)
%
% Description
%
%     dilate is a simple wrapper around imdilate that handles nans.
%
% See also: IMDILATE
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 6. February, 2013




if isa(DEM.Z,'double') || isa(DEM.Z,'single')
    I = isnan(DEM.Z);
    DEM.Z(I) = -inf;
    checknans = true;
else
    checknans = false;
end

DEM.Z = imdilate(DEM.Z,SE);
if checknans
    DEM.Z(I) = nan;
end
