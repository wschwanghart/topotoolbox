function MS = zonalstats(MS,attributes,varargin)

%ZONALSTATS Zonal statistics
%
% Syntax
%
%     MS = zonalstats(MS)
%     MS = zonalstats(MS,{varname1, vargrid1, varfun1, ...
%                         varname2, vargrid1, varfun2, ...})
%     MS = zonalstats(MS,{varname1, vargrid1, varfun1, ...
%                         varname2, vargrid1, varfun2, ...},...
%                        pn,pv,...)
%
% Description
%
%     MS = zonalstats(MS) will add a new field 'area' to the mapping 
%     structure MS that includes the area (calculated by the function 
%     polyarea). 
%
%     MS = zonalstats(MS,{varname1, vargrid1, varfun1, ...
%                         varname2, vargrid2, varfun2, ...}) calculates
%     zonal statistics for GRIDobjs vargrid1, vargrid2, ... using the
%     anonymous functions varfun1, ... and returns the values to a new
%     field in MS entitled varname1, ... . If varfunX is empty, zonalstats
%     will return several fields for the respective vargridX calculated by
%     numerous functions (mean, std, median, percentiles, skewness,
%     kurtosis, etc.).
%
%     MS = zonalstats(...,pn,pv) allows setting a number of options (see
%     below)
%     
%     Note: The function adds a field 'TTID' (TopoToolbox ID) to the output
%     mapping structure.
%
% Input parameters
%
%     MS      polygon mapping structure (e.g. shapefile imported by the
%             function shaperead). The coordinate system must be equal to
%             the grids' projections and should be a projected coordinate
%             system with metric units.
%     {varname1, vargrid1, varfun1, ...} cell array with triplets of
%             variable names, grids and aggregation functions. E.g.
%             {'z' DEM @mean} will create a new field MS.z in MS that 
%             contains the averages of values in DEM within the polygons 
%             defined by MS. If you provide an empty array as varfun1, 
%             zonalstats will calculate several aggregating functions for
%             the respective grid. In this case, zonalstats will create 
%             multiple new fields in MS whose naming starts with varname1
%             followed by underscore _ and abbreviations of the aggregating
%             functions (e.g. z_mean, z_std, z_min, etc.).
%
% Parameter name/value functions
%
%      'overlapping' true or {false}. Set to false if polygons in MS are not
%                    overlapping. This will signicantly increase the speed
%                    of the function.
%      'centroid'    true or {false}. True will add the fields xc and yc to 
%                    MS that contain the x and y coordinates of the
%                    centroids, respectively.
%                    Note that centroid is currently calculated based on
%                    the vertex coordinates. Nonuniform spacing between
%                    vertices will bias the centroid.
%      'area'        true or {false}. True will calculate the area of the
%                    polygons in MS in horizontal units of MS.
%      'perimeter'   true or {false}. True will calculate the perimeter of 
%                    the polygons in MS in horizontal units of MS.
%      'bbox'        true or {false}. True will add eight attributes to the
%                    resulting structure array: lucornerx, lucornery,
%                    llcornerx, llcornery, ...
%      'lucorner'    true or {false} calculates the coordinates of the left 
%                    upper corner of the axis-aligned bounding box of the
%                    region. Same can be used for llcorner, rucorner, etc.
%
% Output arguments
%
%      MS            mapping structure with zonal statistics
%
% 
% See also: polygon2GRIDobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. September, 2020

if nargin == 1
    a = cellfun(@(x,y) getarea(x,y),{MS.X},{MS.Y},'UniformOutput',false);
    [MS.area] = a{:};
    return
end

% check attributes
if mod(numel(attributes),3) ~= 0
    error('additional arguments must come in triplets')
end

% check the attributes
attname      = attributes(1:3:end);
attgrid      = attributes(2:3:end);
attfun       = attributes(3:3:end);
nrattributes = numel(attributes)/3;

p = inputParser;
addParameter(p,'overlapping',false)
addParameter(p,'centroid',false)
addParameter(p,'area',false)
addParameter(p,'perimeter',false)
addParameter(p,'bbox',false);
addParameter(p,'lucorner',false);
addParameter(p,'llcorner',false);
addParameter(p,'rucorner',false);
addParameter(p,'rlcorner',false);

parse(p,varargin{:});
    
% convert strings to functions, if necessary
for r = 1:nrattributes
    if ischar(attfun{r}) && ~isempty(attfun{r})
        attfun{r} = str2func(attfun{r});
    elseif isempty(attfun{r})
        % do nothing
    end
    
    if r > 2
        validatealignment(attgrid{r},attgrid{1});
    end
    
end

if p.Results.area
    a = cellfun(@(x,y) getarea(x,y),{MS.X},{MS.Y},'UniformOutput',false);
    [MS.area] = a{:};
end
if p.Results.centroid
    [xc,yc] = cellfun(@(x,y) getcentroid(x,y),{MS.X},{MS.Y},'UniformOutput',false);
    [MS.xc] = xc{:};
    [MS.yc] = yc{:};
end
if p.Results.perimeter
    per = cellfun(@(x,y) getperimeter(x,y),{MS.X},{MS.Y},'UniformOutput',false);
    [MS.perimeter] = per{:};
end

if p.Results.lucorner  || p.Results.bbox
    xc = cellfun(@(x) min(x),{MS.X},'UniformOutput',false);
    yc = cellfun(@(y) max(y),{MS.Y},'UniformOutput',false);
    [MS.lucornerx] = xc{:};
    [MS.lucornery] = yc{:};
end

if p.Results.rucorner  || p.Results.bbox
    xc = cellfun(@(x) max(x),{MS.X},'UniformOutput',false);
    yc = cellfun(@(y) max(y),{MS.Y},'UniformOutput',false);
    [MS.rucornerx] = xc{:};
    [MS.rucornery] = yc{:};
end
    
if p.Results.llcorner  || p.Results.bbox
    xc = cellfun(@(x) min(x),{MS.X},'UniformOutput',false);
    yc = cellfun(@(y) min(y),{MS.Y},'UniformOutput',false);
    [MS.llcornerx] = xc{:};
    [MS.llcornery] = yc{:};
end

if p.Results.rlcorner  || p.Results.bbox
    xc = cellfun(@(x) max(x),{MS.X},'UniformOutput',false);
    yc = cellfun(@(y) min(y),{MS.Y},'UniformOutput',false);
    [MS.rlcornerx] = xc{:};
    [MS.rlcornery] = yc{:};
end

if ~p.Results.overlapping
    TTID = num2cell(uint32(1:numel(MS)));
    [MS.TTID] = TTID{:};
    label = polygon2GRIDobj(attributes{2},MS,'TTID');
    TTID = num2cell((1:numel(MS)));
    [MS.TTID] = TTID{:};
    stats = regionprops(label.Z,'PixelIdxList');
end



nr = numel(MS);

if nr > 2;
h = waitbar(0,['0 processed, ' num2str(nr) ' remaining']);
wb = true;
else
wb = false;
end
for r = 1:nr;
    
    if p.Results.overlapping
        II = polygon2GRIDobj(attgrid{1},MS(r));
        I  = II.Z;
    else
        I = stats(r).PixelIdxList;
    end
    
    for r2 = 1:nrattributes
        if ~isempty(attfun{r2})
            MS(r).(attname{r2}) = attfun{r2}(double(attgrid{r2}.Z(I)));  
        else
            z = attgrid{r2}.Z(I);
            z = double(z(:));
            z(isnan(z)) = [];
            MS(r).([attname{r2} '_mean']) = mean(z);
            MS(r).([attname{r2} '_std']) = std(z);
            minz = min(z);
            maxz = max(z);
            if isempty(minz); minz = nan; end
            if isempty(maxz); maxz = nan; end;
            MS(r).([attname{r2} '_min']) = minz;
            MS(r).([attname{r2} '_max']) = maxz;
            MS(r).([attname{r2} '_median']) = median(z);
            prc = quantile(z,[.01 .05 .33 .66 .95 .99]);
            MS(r).([attname{r2} '_prc01']) = prc(1);
            MS(r).([attname{r2} '_prc05']) = prc(2);
            MS(r).([attname{r2} '_prc33']) = prc(3);
            MS(r).([attname{r2} '_prc66']) = prc(4);
            MS(r).([attname{r2} '_prc95']) = prc(5);
            MS(r).([attname{r2} '_prc99']) = prc(6);
            MS(r).([attname{r2} '_skew']) = skewness(z);
            MS(r).([attname{r2} '_kurt']) = kurtosis(z);
        end
    end
    if wb
    waitbar(r/nr,h,[num2str(r) ' processed, ' num2str(nr-r) ' remaining']);
    end
end

if wb
close(h)
end




end

function a = getarea(x,y)
I = ~isnan(x);
a = double(polyarea(x(I),y(I)));
end

function [xc,yc] = getcentroid(x,y)

I = ~isnan(x);
[xc,yc] = centroid(polyshape(x(I),y(I)));

% xc = double(mean(x(I)));
% yc = double(mean(y(I)));
end

function p = getperimeter(x,y)
I = ~isnan(x);
xd = diff(x(I));
yd = diff(y(I));

p  = double(sum(hypot(xd,yd)));

end