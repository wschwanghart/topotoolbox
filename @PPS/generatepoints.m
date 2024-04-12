function P = generatepoints(P,varargin)
%GENERATEPOINTS Generate non-random points on stream network
%
% Syntax
%
%     P2 = generatepoints(P)
%     P2 = generatepoints(P,'n',n)
%
% Description
%
%     generatepoints computes a set of points distributed on the
%     streamnetwork stored in P.
%
%     P2 = generatepoints(P) generates points that are equally spaced along
%     20 distance values from the outlet.
%
%     P2 = generatepoints(P,'n',n) generates points that are equally spaced
%     along n distance values from the outlet. Note that P2 will contain
%     more than n points if the underlying stream network is a branched
%     network and/or contains several drainage basins.
%
%     P2 = generatepoints(P,'n',n,'distance',d) generates points that are
%     equally spaced in a user-defined distance d. d must be a node
%     attribute list and strictly monotonically increasing in upstream
%     direction.
%
%     P2 = generatepoints(P,'n',n,'distance',d,'rmoutletpoint',true)
%
%     P2 = generatepoints(P,'type',type) generates points that are either
%     'channelheads', 'outlets', 'confluences', 'bconfluences'. See
%     STREAMobj/streampoi for details. In addition, type can be 'centroid' 
%     which calculates the network centroid.
%
%
% Input arguments
%
%     P     instance of PPS
%     
%     Parameter name/value pairs
%
%     'add'        {false} or true. If true, than points will be added to
%                  the existing point pattern in P.
%
%     For equally spaced points
%
%     'n'          positive integer scalar. Number of distance values. The 
%                  default is 20.
%     'distance'   node attribute list of distance values. The values must
%                  be monotonically increasing in upstream direction. 
%
%     For points of interest
%
%     'type'       see STREAMobj/streampoi for types
%
% Output arguments
%
%     P     instance of PPS
%
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,'minarea',500);
%     P   = PPS(S,'rpois',0.0001);
%     P2  = generatepoints(P,'n',40);
%     plot(P2)   
%
% Example 2
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     A   = flowacc(FD);
%     S   = STREAMobj(FD,'minarea',500);
%     S   = klargestconncomps(S);
%     P   = PPS(S,'rpois',0.0001,'z',DEM);
%     c   = chitransform(S,A);
%     P   = generatepoints(P,'n',10,'distance',c);
%     subplot(1,2,1)
%     plotdz(P)
%     subplot(1,2,2)
%     plotdz(P,'distance',c)
%
%     figure
%     rhohat(P,'covariate',c,'bandwidth',100,'name','\chi')
%
% 
% See also: PPS, PPS/simulate, PPS/random 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2020


p = inputParser;
addParameter(p,'distance',P.S.distance,@(x) isnal(P.S,x));
addParameter(p,'n',20,@(x) validateattributes(x,{'numeric'},{'>',1}));
addParameter(p,'add',false);
addParameter(p,'rmoutletpoint',true);
addParameter(p,'type','');
parse(p,varargin{:});


if isempty(p.Results.type)
    d = p.Results.distance;
    % check if supplied distance is monotonically increasing in upstream direction
    tf = any(gradient(P.S,d) < 0);
    if tf
        error('Distance must be strictly monotonically increasing in upstream direction')
    end
    
    n = p.Results.n;
    
    mind = min(d);
    maxd = max(d);
    
    dvals = linspace(mind,maxd,n+2);
    dvals = dvals(2:end-1);
    
    ix = P.S.ix';
    ixc = P.S.ixc';
    
    I = sparse(d(ixc)) <= sparse(dvals) & sparse(d(ix)) >= sparse(dvals);
    [r,c] = find(I);
    d_to_node = [d(ix(r))-dvals(c)' dvals(c)'-d(ixc(r))];
    [~,which_neighbor] = min(d_to_node,[],2);
    I = which_neighbor == 1;
    
    pp = ix(r(I))';
    pp = [pp; ixc(r(~I))'];
    if p.Results.add
        P.PP = [P.PP; pp];
    else
        P.PP = pp;
    end
    
    if p.Results.rmoutletpoint
        P2 = P;
        P2.PP = [];
        P2 = generatepoints(P2,'type','outlet');
        P.PP(ismember(P.PP,P2.PP)) = [];
    end
    
    
else
    
    if ~iscell(p.Results.type)
        switch lower(p.Results.type)

            case 'centroid'
            P  = generatepoints(P,'type',{'outlets','channelheads'});
            cc = conncomps(P.S);
            cc = getmarks(P,cc);
            P2  = aggregate(P,cc);
            
            if p.Results.add
                P.PP = [P.PP; P2.PP];
            else
                P = P2;
            end
            otherwise
                nal = streampoi(P.S,p.Results.type,'logical');
                P.PP = find(nal);
                
        end
    else
        I = streampoi(P.S,p.Results.type,'logical');
        if p.Results.add
            P.PP = [P.PP; find(I)];
        else
            P.PP = find(I);
        end
    end
end