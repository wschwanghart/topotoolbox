function S = union(varargin)

%UNION merge different instances of STREAMobj into a new instance
%
% Syntax
%
%     S = union(S1,S2,...)
%     S = union(S1,S2,...,FD)
%
% Description
%
%     union combines different instances of STREAMobj into a new STREAMobj.
%
%     union(S1,S2,...) combines all instances into a new STREAMobj. The
%     networks might be subset of each others. All networks must have been
%     derived from the same flow direction object (FLOWobj). If this is not the case, the function might have an
%     unexpected behavior. 
%
%     union(S1,S2,...,FD) takes in addition an instance of FLOWobj as last input
%     variable and it is assumed that all STREAMobjs have been derived
%     from this FLOWobj (e.g., all STREAMobjs are subgraphs of the
%     FLOWobj). This syntax reestablishes the connectivity between adjacent
%     nodes that are connected in the FLOWobj but not in the STREAMobjs.
%
% Input arguments
%
%     S1, S2, ...    several instances of STREAMobj
%     FD             instance of FLOWobj
%
% Output arguments
%
%     S              STREAMobj
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     A  = flowacc(FD);
%     S  = STREAMobj(FD,A>1000);
%     CS = STREAMobj2cell(split(S));
%     S2 = union(CS{:});
%     ax = subplot(2,1,1);
%     plot(S2)
%     S3 = union(CS{:},FD);
%     ax(2) = subplot(2,1,2);
%     plot(S3);
%     linkaxes(ax,'xy')
%
% See also: STREAMobj, FLOWobj, STREAMobj/modify, STREAMobj/trunk
%           STREAMobj/STREAMobj2cell, STREAMobj/intersect
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. October, 2013

narginchk(2,inf)

% are instances spatially aligned?
for r = 1:numel(varargin);
    validatealignment(varargin{1},varargin{r})
end

% is the last input argument a FLOWobj
if isa(varargin{end},'FLOWobj')
    
    % empty GRIDobj
    G = GRIDobj([]);
    % find common properties of F and G and from F to G
    pg = properties(G);
    pf = properties(varargin{end});
    p  = intersect(pg,pf);
    for r = 1:numel(p);
        G.(p{r}) = varargin{end}.(p{r});
    end
    G.Z = false(G.size);
    
    % create a stream raster
    for r = 1:(numel(varargin)-1);
        G.Z(varargin{r}.IXgrid) = true;
    end
    
    % use stream raster to create a new STREAMobj
    S = STREAMobj(varargin{end},G);
    
    
else
    
    
    % vertical concatenate indices
    IXgrid = [];
    ix     = [];
    ixc    = [];
    nnodes = 0;
    for r = 1:numel(varargin)
        IXgrid = [IXgrid;varargin{r}.IXgrid];
        ix     = [ix;  varargin{r}.ix  + nnodes];
        ixc    = [ixc; varargin{r}.ixc + nnodes];
        nnodes = nnodes + numel(varargin{r}.IXgrid);
    end
    
    [IXgrid,~,ic] = unique(IXgrid,'stable');
    
    IX = (1:numel(IXgrid))';
    IX = IX(ic);
    ix = IX(ix);
    ixc = IX(ixc);
    
    ixixc = unique([ix(:) ixc(:)],'rows','stable');
    ix    = ixixc(:,1);
    ixc   = ixixc(:,2);
    
    [ix,ixc] = updatetoposort(nnodes,ix,ixc);
    
    S        = varargin{r};
    S.ix     = ix;
    S.ixc    = ixc;
    S.IXgrid = IXgrid;
    
    [r,c] = ind2sub(S.size,S.IXgrid);
    xy    = double([r c ones(numel(S.IXgrid),1)])*S.refmat;
    S.x = xy(:,1);
    S.y = xy(:,2);
    
end
end 

function [ix,ixc] = updatetoposort(nnodes,ix,ixc)

    nr= nnodes;
    D = digraph(ix,ixc);
    p = toposort(D);
    
    IX = zeros(nr,1,'uint32');
    IX(ix) = ix;
    IXC = zeros(nr,1,'uint32');
    IXC(ix) = ixc;
    
    p   = p(:);
    IX  = IX(p);
    IXC = IXC(p);
    IX  = nonzeros(IX);
    IXC = nonzeros(IXC);
    ix  = double(IX(:));
    ixc = double(IXC(:));
end


