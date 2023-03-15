function [CS,locS,order] = STREAMobj2cell(S,ref,n)

%STREAMOBJ2CELL convert instance of STREAMobj to cell array of stream objects
%
% Syntax
%
%     CS = STREAMobj2cell(S)
%     CS = STREAMobj2cell(S,'outlets')
%     CS = STREAMobj2cell(S,'tributaries')
%     CS = STREAMobj2cell(S,'channelheads')
%     CS = STREAMobj2cell(S,'segments')
%     CS = STREAMobj2cell(S,'outlets',n)
%     CS = STREAMobj2cell(S,'segments',seglength)
%     CS = STREAMobj2cell(S,'labels',labels)
%     CS = STREAMobj2cell(S,'split',ix)
%     CS = STREAMobj2cell(S,'split',P)
%     [CS,locS] = ...
%     [CS,locS,order] = STREAMobj2cell(S,'tributaries');
%
% Description
%
%     STREAMobj2cell divides a STREAMobj into a number of STREAMobj stored
%     in a cell array. 
%
%     STREAMobj2cell(S,'outlets') divides a STREAMobj into its strongly
%     connected components. This means that individual STREAMobjs are
%     derived as individual trees of the stream network, i.e. individual
%     drainage basins. This is the default. In this case, CS has as many
%     elements as there are outlets in the stream network. 
%
%     STREAMobj2cell(S,'outlets',n) takes the n largest connected 
%     components and stores these in the cell array CS. In this case, CS
%     has n elements.
%     
%     STREAMobj2cell(S,'channelheads') derives individual STREAMobj as
%     single streams emanating from each channelhead. In this case, CS has
%     as many elements as there are channelheads in the stream network.
%
%     STREAMobj2cell(S,'tributaries') derives individual STREAMobj as
%     tributaries. A stream is a tributary until it reaches a stream with a 
%     longer downstream flow distance.   
%
%     STREAMobj2cell(S,'segments') splits the stream network into
%     individual reaches. By default, the maximum reach length is 20 the
%     cellsize. A third arguments controls the segment length. Note that
%     the network is always split at river junctions.
%
%     STREAMobj2cell(S,'labels',labels) takes a node-attribute list with
%     labels as third input argument. Nodes with common labels will be
%     placed into the same element in CS.
%
%     STREAMobj2cell(S,'split',ix) splits the stream network at cells with
%     the linear index ix. 
%
%     STREAMobj2cell(S,'split',P) takes an instance of PPS and splits the
%     stream network S at points stored in P. Ideally, the underlying
%     stream network in P should be S.
%     
% Input arguments
%
%     S          instance of STREAMobj
%     ref        reference for deriving individual STREAMobj. Either 
%                'outlets' (default) or 'channelheads', or 'tributaries', 
%                or 'segments'.
%     n          if ref is 'outlets', n determines the number of n largest 
%                trees to be placed in elements of CS.
%     seglength  maximum segment length if ref is 'segments'
%     ix         linear index at which the stream network is split. The
%                cells with the linear index ix must be part of the stream
%                network
%
% Output arguments
%
%     CS    cell array with elements of CS being instances of STREAMobj
%     locS  cell array with linear indices into node attribute lists of S
%     order output argument only if ref is set to 'tributaries'. Vector
%           with numel(CS) elements where each element refers to an order 
%           of the tributaries. 
%
% Example 1: Place stream networks of individual basins into a cell array
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','carve');
%     S   = STREAMobj(FD,'minarea',1000);
%     z   = getnal(S,DEM);
%     [CS,locS] = STREAMobj2cell(S);
%     plotdz(CS{21},z(locS{21}))
%
% Example 2: Create random river segments and place them in a cell array
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S,1);
%     P = PPS(S,'rpois',0.0001,'z',DEM);
%     plot(P)
%     CS = STREAMobj2cell(S,'split',P);
%     figure;
%     hold on; cellfun(@plot,CS)
%
% See also: FLOWobj2cell, STREAMobj/split, PPS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. March, 2022

if nargin == 1
    ref = 'outlets';
    getall = true;
    n   = inf;
elseif nargin == 2
    ref = validatestring(ref,{'outlets','channelheads','tributaries','segments'},'STREAMobj2cell','ref',2);
    getall = true;
    n   = inf;
    seglength = 20.*S.cellsize;
elseif nargin == 3
    ref = validatestring(ref,{'outlets','segments','labels','split'},'STREAMobj2cell','ref',2);
    switch ref
        case {'outlets','segments'}
            validateattributes(n,{'numeric'},{'>',1},'STREAMobj2cell','n',3);
        case 'labels'
            if ~isnal(S,n)
                error('Array must be a node attribute list')
            end
    end
    getall = false;
    seglength = n;
end

switch lower(ref)
    case 'outlets'
        
        nrc = numel(S.x);
        M = sparse(double(S.ix),double(S.ixc),true,nrc,nrc);
        
        [~,p,~,r] = dmperm(M | M' | speye(nrc));
               
        nc = length(r) - 1;
        if getall || nc < n
            % label matrix
            L = zeros(nrc,1);
            for tt = 1:nc
                L(p(r(tt):r(tt+1)-1)) = tt;
            end
        else
            nc = n;
            L = zeros(nrc,1)+nc;
            
            [~,dd] = sort(diff(r),'descend');
            
            counter = 1;

            for tt = dd(1:nc) %1:min(nc,k);
                L(p(r(tt):r(tt+1)-1)) = counter;
                counter = counter + 1;
            end
            
        end
        
        % put each individual tree into an own element in CS.
        CS = cell(1,nc);
        
        % adapt new STREAMobj to the reduced network
        LL = L;
        
        if nargout == 2
            locS = cell(1,nc);
        end
        
        for r = 1:nc
            CS{r} = S;
            L     = LL==r;
            I     = L(CS{r}.ix);
            CS{r}.ix  = CS{r}.ix(I);
            CS{r}.ixc = CS{r}.ixc(I);
            
            IX    = cumsum(L);
            CS{r}.ix  = IX(CS{r}.ix);
            CS{r}.ixc = IX(CS{r}.ixc);
            
            CS{r}.x   = CS{r}.x(L);
            CS{r}.y   = CS{r}.y(L);
            CS{r}.IXgrid   = CS{r}.IXgrid(L);
            
            if nargout == 2
                locS{r} = find(L);
            end
                
            
        end
        
        
        
    case 'channelheads'
        ch     = streampoi(S,'channelheads','logical');
        ixcix  = zeros(numel(S.IXgrid),1);
        ixcix(S.ix) = 1:numel(S.ix);
        
        ixchannel = find(ch);
        nc = numel(ixchannel);
        
        CS = cell(1,nc);
        if nargout == 2
            locS = cell(1,nc);
        end
        for r = 1:nc
            IX = ixchannel(r);
        
            c = 1;
            while ixcix(IX) ~= 0
                c = c+1;
                IX(c) = S.ixc(ixcix(IX(end)));
            end
            
            L = false(size(S.IXgrid));
            L(IX) = true;
            
            % adapt new STREAMobj to the reduced network
            
            CS{r} = S;
            I     = L(CS{r}.ix);
            CS{r}.ix  = CS{r}.ix(I);
            CS{r}.ixc = CS{r}.ixc(I);
            
            IX    = cumsum(L);
            CS{r}.ix  = IX(CS{r}.ix);
            CS{r}.ixc = IX(CS{r}.ixc);
            
            CS{r}.x   = CS{r}.x(L);
            CS{r}.y   = CS{r}.y(L);
            CS{r}.IXgrid   = CS{r}.IXgrid(L);
            
            if nargout == 2
                locS{r} = find(L);
            end
            
            
            
        end
        
    case 'tributaries'
        if nargout < 3
            CS = tributaries(S);
        else
            CS = tributariesinclorder(S);
            order = CS(2:2:end);
            order = horzcat(order{:});
            CS = CS(1:2:end);
        end
        
        if nargout >= 2
            locS = cell(size(CS));
            for r = 1:numel(CS)
                [~,locS{r}] = ismember(CS{r}.IXgrid,S.IXgrid);
            end
        end
        
        return
        
    case {'segments','labels'}
        
        switch ref
            case 'segments'
                lab = labelreach(S,'seglength',seglength);
            case 'labels'
                [~,~,lab] = unique(seglength);
        end
        nc  = max(lab);
        CS  = cell(1,nc);
        if nargout == 1
            parfor r = 1:numel(CS)
                CS{r} = subgraphix(S,lab==r);
            end
        else
            locS = cell(1,nc);
            parfor r = 1:numel(CS)
                [CS{r},locS{r}] = subgraphix(S,lab==r);
            end
        end
        
    case 'split'
        
        if isa(n,'PPS')
            n = n.S.IXgrid(n.PP);
        end
            
        
        ix = n;
        ix = unique(ix);
        [I,locb]  = ismember(S.IXgrid,ix);
        
        % Remove outlet pixels if they are set to true
        outlets = streampoi(S,'outlet','logical');
        I(outlets) = false;
        
        I2 = false(size(I));
        I2(S.ix) = I(S.ixc);
        
        I  = I2;

        I(outlets) = true;
        I(streampoi(S,'chan','logical')) = false;
        
        label = zeros(size(S.x));
        label(I) = 1:nnz(I);
        
        for r = numel(S.ix):-1:1
            if label(S.ix(r)) == 0
                label(S.ix(r)) = label(S.ixc(r));
            end
        end
        
        if nargout == 1
            CS = STREAMobj2cell(S,'label',label);
        else
            [CS,locS] = STREAMobj2cell(S,'label',label);
        end
        return  
end

% check for validity of Ss
valid = true(1,nc);
for r = 1:nc
    if numel(CS{r}.x) == 1
        valid(r) = false;
    end
end
CS = CS(valid);

if nargout == 2
    locS = locS(valid);
end
end


%% Recursively scan for tributaries
function Ctribs = tributaries(S)
% Recursive extraction of tributaries


C = STREAMobj2cell(S);
Ctribs = cell(0);

for r = 1:numel(C)
    St = trunk(C{r});
    S2 = modify(C{r},'tributaryto2',St);
    if isempty(S2.ix)
        % do nothing
        Ctribs   = [Ctribs {St}];
    else
        Ctribs = [Ctribs {St} tributaries(S2)];
    end
end
end



function Ctribs = tributariesinclorder(S,orderin)
% Recursive extraction of tributaries with order


C = STREAMobj2cell(S);
Ctribs = cell(0);

if nargin == 1
    orderin = 1;
end
    

for r = 1:numel(C)
    St = trunk(C{r});
    S2 = modify(C{r},'tributaryto2',St);
    if isempty(S2.ix)
        % do nothing
        Ctribs   = [Ctribs {St orderin}];
    else
        Ctribs = [Ctribs {St orderin} tributariesinclorder(S2, orderin+1)];
    end
end
end

function [S,locb] = subgraphix(S,nal)


if all(nal)
    % do nothing
    return
end

if nargout == 2
    IXgrid_old = S.IXgrid;
end

I = nal(S.ix);

S.ix  = S.ix(I);
S.ixc = S.ixc(I);

% nal  = nal;
nal(S.ixc) = true;
IX    = cumsum(nal);

S.ix  = IX(S.ix);
S.ixc = IX(S.ixc);

S.x   = S.x(nal);
S.y   = S.y(nal);
S.IXgrid   = S.IXgrid(nal);

if nargout == 2
    [~,locb] = ismember(S.IXgrid,IXgrid_old);
end

end

