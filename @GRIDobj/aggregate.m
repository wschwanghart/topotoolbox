function C = aggregate(DEMlowres,DEMhighres,aggfun)

%AGGREGATE resampling a GRIDobj using aggregation
%
% Syntax
%
%     C = aggregate(B,A)
%     C = aggregate(B,xy)
%     C = aggregate(B,xyz)
%     C = aggregate(...,aggfun)
%
% Description
%
%     This function resamples the grid A to C to match the extent and
%     resolution of grid B. B must spatially overlap with A and must have a
%     coarser resolution. By default, aggregate uses the mean to calculate
%     new values in grid C, but aggregate takes any other function that takes 
%     a vector and returns a scalar (e.g. median, std, ...).
%
%     Values to be aggregated can also be supplied as list of coordinates
%     (and attributes). This is particularly useful if point density for
%     each pixel is greater than one. 
%
% Input arguments
%
%     A         high resolution GRIDobj
%     B         low resolution GRIDobj. A and B must have the same 
%               coordinate system
%     xy        two column matrix with coordinates. aggfun is @numel.
%     xyz       three column matrix with coordinates and third column being
%               some attribute (e.g. elevation, ...)
%     aggfun    anonymous function that aggregate values. The function must
%               take a vector and return a scalar. The default is @mean.
%               Other possible functions are @numel to obtain counts,
%               @median, @std, ...
%
% Output arguments
%
%     C         GRIDobj with same extent and resolution as B and aggregated
%               values derived from A
%
%
% See also: GRIDobj/resample, accumarray
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 12. July, 2018

if nargin == 2
    aggfun = @mean;
end

if isa(DEMhighres,'GRIDobj')
    % A is a high res GRIDobj
    [x,y] = getcoordinates(DEMhighres);
    IX    = zeros(DEMhighres.size);
    [X,Y] = getcoordinates(DEMlowres);

    sy = size(y);
    for r = 1:numel(x)
        IX(:,r) = coord2ind(X,Y,repmat(x(r),sy),y);
    end
    
    Z = DEMhighres.Z(:);
    IX = IX(:);
    inan = isnan(IX);
    IX(inan) = [];
    Z(inan) = [];
    
else
    % A is a list of coordinates and attributes
    xyz = DEMhighres;
    IX  = coord2ind(DEMlowres,xyz(:,1),xyz(:,2));
    inan = isnan(IX);
    xyz(inan,:) = [];
    IX(inan) = [];
    
    if size(xyz,2) == 3
        Z = xyz(:,end);
    else
        Z = true; %true(size(IX));
        aggfun = @numel;
    end
end
    
if isempty(Z)
    C = GRIDobj(DEMlowres);
    return
end
    
z = accumarray(IX,Z,[prod(DEMlowres.size) 1],aggfun);
C = DEMlowres;
C.Z = reshape(z,DEMlowres.size);
% C.Z(C.Z == 0) = nan;


