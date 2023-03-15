function [v,nIX,CP] = extractvaluesaroundpoints(P,z,varargin)

%EXTRACTVALUESAROUNDPOINTS Extract values around points 
%
% Syntax
%
%     v = extractvaluesaroundpoints(P,z)
%     v = extractvaluesaroundpoints(p,z,'pn',pv,...)
%     [v,nIX,CP] = ...
%
% Description
%
%     extractvaluesaroundpoints extracts values from a GRIDobj or
%     node-attribute list z along the stream network around the points in
%     P. P is an instance of PPS and stores the stream network and points
%     located on the stream network. 
%
%     The function can be, for example, used to extract river gradients (or
%     ksn values) upstream and downstream of knickpoints, the difference of
%     which can be used as a metric of knickpoint prominence.
%
% Input arguments
%
%     P     instance of PPS
%     z     GRIDobj or node-attribute list
%     
%     Parameter name/value pairs
%
%     'direction'   {'both'}, 'upstream', 'downstream' or 'bothundirected'
%                   'both' extracts the stream network in upstream and
%                   downstream direction. It won't include downstream
%                   tributaries, however. 'bothundirected' will include all
%                   stream segments, irrespective of their direction.
%     'dfrompoint'  distance from points. If 'direction' is 'both', then
%                   'dfrompoint' accepts a two-element vector where there
%                   first element refers to the distance downstream from
%                   the point, and the second element refers to the
%                   distance upstream from the point.
%     'aggfun'      anonymous aggregation function which takes a vector and 
%                   returns a scalar {@mean}. 
%     'distance'    By default, distances between are calculated as the
%                   euclidean distance between stream nodes. Yet, other
%                   distance metrics (e.g. chi) can be used as well. Note
%                   that the value in dfrompoint must be given in the same
%                   distance units as those provided by here.
%    
% Output arguments
%
%     v       extracted value
%     nIX     cell array of linear indices into the grid.
%     CP      cell array of PPS objects, one for each point, with the river
%             network extracted. 
%    
% Example: Extract ksn values upstream and downstream of knickpoints
%
%     DEM  = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S,1);
%     [~,kp] = knickpointfinder(S,DEM,'tol',30,'verbose',false,...
%         'plot',false);
%     P = PPS(S,'PP',kp.IXgrid,'z',DEM);
%     ks = ksn(S,DEM,flowacc(FD));
%     
%     ksd = extractvaluesaroundpoints(P,ks,'direction','down',...
%         'dfrompoint',200);
%     ksu = extractvaluesaroundpoints(P,ks,'direction','up',...
%         'dfrompoint',200);
%     scatter(ksd,ksu)
%     xlabel('Ksn downstream')
%     ylabel('Ksn upstream')
%     refline(1,0)
%     box on
%
% See also: PPS, PPS/getmarks, nearest
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 20. September, 2022

p = inputParser;
p.FunctionName = 'extractvaluesaroundpoints';

addParameter(p,'direction','both',@(x) ischar(validatestring(x,{'both','upstream','downstream','bothundirected'})))
addParameter(p,'dfrompoint',10*P.S.cellsize)
addParameter(p,'distance',P.S.distance)
addParameter(p,'aggfun',@mean)

parse(p,varargin{:})

if isa(z,'GRIDobj')
    z = double(getnal(P.S,z));
else
    if isnal(P.S,z)
        z = double(z);
    else
        error('Cannot handle second input. It must be a GRIDobj or node-attribute list.')
    end
end

% get distance
d = p.Results.distance;
% and calculate inter-node distances as weights
w = abs(d(P.S.ix) - d(P.S.ixc));

direction = validatestring(p.Results.direction,{'both','upstream','downstream'});
switch direction
    case 'downstream'
        G = digraph(P.S.ix, P.S.ixc,w);
    case 'upstream'
        G = digraph(P.S.ixc, P.S.ix, w);
    case 'bothundirected'
        G = graph([P.S.ix; P.S.ixc],[P.S.ixc; P.S.ix], [w;w]);
    case 'both'
        % downstream
        G1 = digraph(P.S.ix, P.S.ixc,w);
        d1 = p.Results.dfrompoint(1);   
        % upstream
        G2 = digraph(P.S.ixc, P.S.ix, w);
        if numel(p.Results.dfrompoint) == 2
            d2 = p.Results.dfrompoint(2);
        else
            d2 = d1;
        end
end

% Calculate nearest neighbor indices using cellfun
pIX = num2cell(P.PP);

if nargout >= 2
    switch direction 
        case 'both'

            nIX = cellfun(@(s) ...
                [nearest(G1,s,d1);...
                s; ...
                nearest(G2,s,d2)],...
                pIX,'UniformOutput',false);


        otherwise
            nIX = cellfun(@(s) [nearest(G,s,p.Results.dfrompoint(1));s],pIX,'UniformOutput',false);
    end
    v   = cellfun(@(ix) p.Results.aggfun(z(ix)),nIX);
    
elseif nargout == 1
    switch direction 
        case 'both'
            v  = cellfun(@(s) p.Results.aggfun(z(...
                [nearest(G1,s,d1);...
                s; ...
                nearest(G2,s,d2)])),...
                pIX,'UniformOutput',true);

        otherwise
            v   = cellfun(@(s) p.Results.aggfun(z([nearest(G,s,p.Results.dfrompoint(1));s])),...
                pIX,'UniformOutput',true);
    end
end

% Extract stream network around points
CP = cell(npoints(P),1);
if nargout >= 3
    for r = 1:npoints(P)
        s = logical(getnal(P.S));
        s(nIX{r}) = true;
        
        if isempty(P.z)
            S = subgraph(P.S,s);
            CP{r} = PPS(S,'PP',P.S.IXgrid(P.PP(r)));
        else
            [S,locb] = subgraph(P.S,s);
            CP{r} = PPS(S,'PP',P.S.IXgrid(P.PP(r)),'z',P.z(locb));
        end
    end
end

if nargout >= 2
    nIX = cellfun(@(ix) P.S.IXgrid(ix),nIX,'UniformOutput',false);
end
