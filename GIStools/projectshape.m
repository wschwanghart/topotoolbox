function [S] = projectshape(S,GRID,varargin)

% PROJECTSHAPE changes the projection of a geographic data structure
%
% Syntax
%
%     MS2 = projectshape(MS,GRID)
%     MS2 = projectshape(MS,GRID,mstruct)
%
% Description
%     
%     projectshape(MS,GRID) transforms the lat,lon of the features in the
%     geographic data structure MS (as obtained from importing a shapefile 
%     using shaperead) to the x,y units of the GRIDobj GRID.
%
%     projectshape(MS,GRID,mstruct) uses the map projection structure
%     mstruct to first convert the x,y units of the features in the
%     geographic data structure MS (as obtained from importing a shapefile 
%     using shaperead) before converting to the x,y units of the GRIDobj 
%     GRID.
% 
%     Note that projectshape requires the Mapping Toolbox
%
% Input arguments
%
%     MS      geographic data structure
%     GRID    GRIDobj, e.g., a DEM
%     mstruct map projection structure
%     
% Output arguments
%
%     MS2     geographic data structure
%
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: 18. September, 2017

[r,c] = size(S(1).X);
if r>c
    X = vertcat(S.X);
    Y = vertcat(S.Y);
else
    X = horzcat(S.X);
    Y = horzcat(S.Y);
end

if min(X)<-180 || max(X)>180 || min(Y)<-90 || max(Y)>90
    
    % Shapefile comes with projected coordinates
    if nargin==3
        mstruct = varargin{1};
        if strcmp(mstruct.mapprojection,'utm')
            % Convert coordinates to lat,lon
            for i = 1 : length(S)
                x = S(i).X;
                y = S(i).Y;
                [lat,lon] = minvtran(mstruct,x,y);
                S(i).X = lon;
                S(i).Y = lat;
            end
        else
            error('Geographic data structure comes with unknown map projection')
        end
    else
        error('No map projection structure found. Conversion of coordinates to lat,lon failed.')
    end
    
end



try
    
    % GRIDobj is projected
    p = projcrs(GRID.georef.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey);
    for i = 1 : length(S)
        lon = S(i).X;
        lat = S(i).Y;
        %[x,y] = mfwdtran(GRID.georef.mstruct,lat,lon);
        [x,y] = projfwd(p,lat,lon);
        S(i).X = x;
        S(i).Y = y;
    end
    
end





