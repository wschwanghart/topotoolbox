function zs = aggregate(S,DEM,varargin)

%AGGREGATE Summarizing values within segments of the stream network
%
% Syntax
%
%     as = aggregate(S,A)
%     as = aggregate(S,a)
%     as = aggregate(...,pn,pv,...)
%
% Description
%
%     Values along stream networks are frequently affected by scatter.
%     This function removes scatter by averaging over values in drainage
%     basins, reaches between confluences, reaches of equal length, or
%     reaches whose endpoints are defined by locations on the stream 
%     network.
%
% Input parameters
%
%     S        STREAMobj
%     A        GRIDobj
%     a        node attribute list
%
% Parameter name/value pairs
%
%     'method'     {'reach'}, 'betweenconfluences','drainagebasins',
%                  'locations'
%                  'reach': values in A are aggregated in river reaches
%                  with segment length as provided in the option
%                  'seglength'
%                  'betweenconfluences': values in A are aggregated in
%                  river sections between confluences.
%                  'drainagebasins': values in A are aggregated in
%                  individual drainagebasins.
%                  'locations': values in A are aggregated between
%                  locations which are provided by the linear index given
%                  in the option 'ix'. See example 2.
%     'split'      {false} or true. True will identify individual drainage
%                  basins and process each individually in parallel 
%                  (requires the parallel processing toolbox).
%     'ix'         if 'method' is 'locations', linear indices of split
%                  locations must be provided. Linear indices must lie on
%                  the stream network (see STREAMobj/snap2stream).
%     'seglength'  segment length (default 10*S.cellsize)
%     'aggfun'     anonymous function as aggregation function. Default is
%                  @mean. The function must take a vector and return a 
%                  scalar (e.g. @max, @min, @std, @(x) prctile(x,25), ...)
%
% Output parameters
%
%     as     node attribute list with aggregated values
%
% Example 1: Calculating ksn values within segments of 2000 m
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     g = gradient(S,imposemin(S,DEM));
%     A = flowacc(FD);
%     a = getnal(S,A)*(DEM.cellsize^2);
%     k = ksn(S,DEM,A);
%     k = aggregate(S,k,'seglength',2000);
%     plotc(S,k)
%     axis image
%     colormap(jet)
%     colorbar
%
% Example 2: Calculating ksn values in river reaches separated by 
%            knickpoints
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     [~,kp] = knickpointfinder(S,DEM,'tol',20,'split',false,...
%                               'plot',false,'verbose',false);
%     A = flowacc(FD);
%     k  = ksn(S,DEM,A);
%     kk = aggregate(S,k,'method','locations','ix',kp.IXgrid);
%     plotc(S,kk)
%     hold on
%     plot(kp.x,kp.y,'+k')
%     colormap(turbo)
%     colorbar
%
% See also: STREAMobj/labelreach, STREAMobj/smooth
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 26. May, 2023


% check and parse inputs
narginchk(2,inf)

p = inputParser;
p.FunctionName = 'STREAMobj/aggregate';
addParameter(p,'method','reach'); 
addParameter(p,'split',false);
addParameter(p,'ix',[]);
% parameters for block and reach
addParameter(p,'seglength',S.cellsize*11);
addParameter(p,'aggfun',@mean);

parse(p,varargin{:});
method = validatestring(p.Results.method,...
    {'betweenconfluences','reach','drainagebasins','locations'});

% get node attribute list with elevation values
if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    z = getnal(S,DEM);
elseif isnal(S,DEM)
    z = DEM;
else
    error('Imcompatible format of second input argument')
end

% run in parallel if wanted
if p.Results.split
    params       = p.Results;
    params.split = false;
    [CS,locS]    = STREAMobj2cell(S);
    Cz           = cellfun(@(ix) z(ix),locS,'UniformOutput',false);
    Czs          = cell(size(CS));
    parfor r = 1:numel(CS)
        Czs{r} = aggregate(CS{r},Cz{r},params);
    end
    
    zs = nan(size(z));
    for r = 1:numel(CS)
        zs(locS{r}) = Czs{r};
    end
    return
end

%% Aggregating starts here

switch method
    case 'betweenconfluences'
        S2 = split(S);
        [c,n] = conncomps(S2);
        za = accumarray(c,z,[n 1],p.Results.aggfun,nan,false);
        zs = za(c);
    case 'reach'
        label = labelreach(S,'seglength',p.Results.seglength);   
        zm    = accumarray(label,z,[max(label) 1],p.Results.aggfun,nan);
        zs    = zm(label);
    case 'drainagebasins'
        label = conncomps(S);
        zm    = accumarray(label,z,[max(label) 1],p.Results.aggfun,nan);
        zs    = zm(label);
    case 'locations'

        S2 = split(S,p.Results.ix);        
        [c,n]  = conncomps(S2);
        za = accumarray(c,z,[n 1],p.Results.aggfun,nan,false);
        zs = za(c);
end

zs = nal2nal(S,S2,zs);

end
        
