function S = intersect(varargin)

%INTERSECT intersect different instances of STREAMobj 
%
% Syntax
%
%     S = intersect(S1,S2,...)
%
% Description
%
%     intersect combines different instances of STREAMobj into a new STREAMobj.
%     All STREAMobjs must have been derived
%     from the same FLOWobj (e.g., all STREAMobjs are subgraphs of the
%     FLOWobj). If this is not the case, the function might have an
%     unexpected behavior. In terms of set theory, the functions produces
%     an intersection of the nodes in all instances of S.
%
% Input arguments
%
%     S1, S2, ...    several instances of STREAMobj
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
%     % get trunk streams of streamorder 3
%     Sn = intersect(trunk(S),modify(S,'streamorder',3));
%     plot(S,'r')
%     hold on
%     plot(Sn,'b')
%     legend('Stream network','Trunk rivers with streamorder = 3')
%     
%
% See also: STREAMobj, FLOWobj, STREAMobj/modify, STREAMobj/trunk
%           STREAMobj/STREAMobj2cell, STREAMobj/union
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. October, 2013


% are instances spatially aligned?
for r = 1:numel(varargin);
    validatealignment(varargin{1},varargin{r})
end

n = numel(varargin);
nrc = prod(varargin{1}.size);
M = sparse(nrc,nrc);
for r = 1:n
    ix  = varargin{r}.IXgrid(varargin{r}.ix);
    ixc = varargin{r}.IXgrid(varargin{r}.ixc);
    M = M + sparse(ix,ixc,ones(size(ix)),nrc,nrc);
end

I = M==n;
[ix,ixc] = find(I);
[IXgrid,~,ic] = unique([ix;ixc]);
ix = ic(1 : numel(ix))';
ixc = ic((numel(ixc)+1) : end)';

[ix,ixc] = updatetoposort(numel(IXgrid),ix,ixc);

S = varargin{1};

S.ix     = ix;
S.ixc    = ixc;
S.IXgrid = IXgrid;

[r,c] = ind2sub(S.size,S.IXgrid);
xy    = double([r c ones(numel(S.IXgrid),1)])*S.refmat;
S.x = xy(:,1);
S.y = xy(:,2);





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

