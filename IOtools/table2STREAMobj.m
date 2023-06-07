function [S,z] = table2STREAMobj(t,DEM,varargin)

%table2STREAMobj Convert point measurements in table to STREAMobj
%
% Syntax
%
%      S = table2STREAMobj(t,DEM)
%      S = table2STREAMobj(t,DEM,'x','xvar','y','yvar')
%      [S,z] = table2STREAMobj(t,DEM,'x','xvar','y','yvar','z','zvar', ...
%                              'direction','up','aggfun',@aggfun)
%
% Description
%
%      The function takes a table which includes horizontal coordinates 'x'
%      and 'y' that describe a path across the grid DEM and returns a
%      STREAMobj S. STREAMobjs are usually retrieved from FLOWobjs and can 
%      be considered as subgraphs of the flow network. Commonly, S returned
%      by this function will not be a subgraph and thus there are
%      limitations to work with it.
%
%      Also note that the function can only handle a single reach, hence a
%      network with only one channelhead.
%
% Input arguments
%
%      t      table
%      DEM    GRIDobj
%      
%      Parameter name/value pairs
%
%      'x'      name of the variable in t which contains the x-coordinates.
%               Default is 'x'.
%      'y'      name of the variable in t which contains the y-coordinates.
%               Default is 'y'.
%      'z'      name of the variable in t which contains the elevation.
%               Default is an empty character array (''). If elevation
%               values are provided, the function will return a node
%               attribute list with elevation values interpolated along the
%               path.
%      'direction' direction, at which the measurements where taken.
%               Default is 'up' and means that measurements were taken 
%               starting from the outlet up to the channelhead. 'down'
%               means that the measurements were taken in the opposite
%               direction.
%      'aggfun' Aggregation function that is applied to elevation values if
%               several  coordinates are associated with the same raster
%               pixel. Default is @mean.
%
% Output arguments
%
%      S      STREAMobj
%      z      node-attribute list with elevation values. If option 'z' is
%             not set, z will be a vector of zeros.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     
%     t = table;
%     t.x = [400000; 390000; 380000; 377000];
%     t.y = [3789000; 3800000; 3789000; 3804000];
%     t.z = [2000; 1500; 1200; 1000];
%     
%     [S,z] = table2STREAMobj(t,DEM,'x','x','y','y','z','z','aggfun',@max);
%     
%     subplot(1,2,1)
%     imageschs(DEM)
%     hold on
%     plot(S,'w')
%     subplot(1,2,2)
%     plotdz(S,z)
%     hold on
%     plotdz(S,DEM)
%
% See also: line2GRIDobj, getgriddedline
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 07. June, 2023



% parse inputs
varnames = t.Properties.VariableNames;
p = inputParser;
addParameter(p,'x','x',@(v) ismember(v,varnames))
addParameter(p,'y','y',@(v) ismember(v,varnames))
addParameter(p,'z','',@(v) ismember(v,varnames) || isempty(v))
addParameter(p,'direction','up',@(x) ischar(validatestring(x,{'up','down'})))
addParameter(p,'aggfun','first')
parse(p,varargin{:})


% We go in downstream direction. If option is up, then we just flip the
% table
switch lower(p.Results.direction)
    case 'up'
        t = flipud(t);
end

% retrieve data from table
x = t.(p.Results.x);
y = t.(p.Results.y);
if ~isempty(p.Results.z)
    z = t.(p.Results.z);
else
    z = zeros(size(x));
end

% map to grid
ix = coord2ind(DEM,x,y);

% Are all ix unique. If not, an aggregation function is used to map
% duplicate values to cells.

[ixun,ia,ib] = unique(ix,'stable');
if ischar(p.Results.aggfun) || isstring(p.Results.aggfun)
switch lower(p.Results.aggfun)
    case 'first'
        
        x = x(ia);
        y = y(ia);
        z = z(ia);
    otherwise
        error('first or function handle is currently the only option.')
end

else
        if isa(p.Results.aggfun,'function_handle')
            x = x(ia);
            y = y(ia);
            z = accumarray(ib,z,[numel(ixun) 1],p.Results.aggfun);
        end
end

% Get path on grid (take from line2GRIDobj)
IX = getgriddedline(DEM,x,y);
[x,y] = ind2coord(DEM,IX);

S = STREAMobj;
S.IXgrid = IX;
S.x = x;
S.y = y;
S.ix = (1:(numel(x)-1))';
S.ixc = (2:(numel(x)))';

S.size = DEM.size;
S.refmat = DEM.refmat;
S.georef = DEM.georef;
S.cellsize = DEM.cellsize;

if nargout == 2
    zp = z;
    z = getnal(S)*nan;
    [~,ix] = ismember(ixun,IX);
    z(ix) = zp;
    z = inpaintnans(S,z);
    % z = interp(S,S.distance(~isnan(z)),z(~isnan(z)));
end


