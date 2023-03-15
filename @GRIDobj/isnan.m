function I = isnan(DEM)

% returns array elements that are NaNs as logical grid
%
% Syntax
%
%     I = isnan(DEM)
%
% Description
%
%     overloaded isnan for GRIDobj. 
%
% See also: isnan
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 20. February, 2013

I = DEM;
I.Z = isnan(DEM.Z);
I.name = [DEM.name ' (isnan)'];