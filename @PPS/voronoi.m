function [nearestnode,d] = voronoi(P,varargin)

%VORONOI nearest neighbor search on a stream network
%
% Syntax
%
%     [v,d] = voronoi(P)
%
% Description
%
%     voronoi returns a node-attribute list v with labels so that each node
%     with the label i is closest to the point P.PP(i). The node-attribute 
%     list d contains the distance of each node to its closest point.
%
% Input arguments
%
%     P      Instance of PPS
%
%     Parameter name/value pairs
%
%     'd3d'         calculate distances in 3d. Requires that P has a 
%                   z-property
%     'distance'    define custom node-attribute list with distance values
%     'shuffle'     shuffle labels randomly. {false} or true. If shuffled,
%                   there is no association to the points in P anymore.
%
% Output arguments
%
%     v      node-attribute list with labels
%     d      node-attribute list with distances
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     S = removeshortstreams(S,100);
%     S = clean(S);
%     P = PPS(S,'rpois',0.001,'z',DEM);
%     [v,d] = voronoi(P);
%     plotc(P.S,v)
%     hold on
%     plotpoints(P,'marks',1:npoints(P))
%     hold off
%     colormap(lines)
% 
% See also: PPS, PPS/npoints 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019



% Check input arguments
p = inputParser;
p.FunctionName = 'PPS/voronoi';
addParameter(p,'distance',[]);
addParameter(p,'d3d',false);
addParameter(p,'shuffle',true);
% Parse
parse(p,varargin{:});


ix  = P.S.ix;
ixc = P.S.ixc;

nodes = zeros(size(P.S.x));
nodes(P.PP) = 1:npoints(P);

% get inter-node distances
dedge = nodedistance(P,'d3d',p.Results.d3d,'val',p.Results.distance);
    
nearestnode = nodes;
d           = inf(size(P.S.x));
d(P.PP)     = 0;

% Find nearest node and calculate distance in downstream direction
for r = 1:numel(ix)
    if nearestnode(ix(r))  && ~nearestnode(ixc(r))
        nearestnode(ixc(r)) = nearestnode(ix(r));
        d(ixc(r)) = min(d(ix(r)) + dedge(r),d(ixc(r)));
    elseif nearestnode(ix(r)) && nearestnode(ixc(r))
        dtest = d(ix(r)) + dedge(r);
        if d(ixc(r)) > dtest
            d(ixc(r)) = dtest;
            nearestnode(ixc(r)) = nearestnode(ix(r));
        end
    end
end

% Find nearest node and calculate distance in upstream direction
dup = inf(size(P.S.x));
dup(P.PP) = 0;
for r = numel(ix):-1:1
    
    dup(ix(r)) = dup(ixc(r)) + dedge(r);
    
    if dup(ix(r)) <= d(ix(r))
        d(ix(r)) = dup(ix(r));
        nearestnode(ix(r)) = nearestnode(ixc(r));
        
    else %dup(ix(r)) > d(ix(r))
        dup(ix(r)) = d(ix(r));
    end
        
end        

if p.Results.shuffle
    [uniqueL,~,ix] = unique(nearestnode);
    uniqueLS = randperm(numel(uniqueL));
    nearestnode = uniqueLS(ix);
    nearestnode = nearestnode(:);
end
    
