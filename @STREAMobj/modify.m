function [Sout,nalix] = modify(S,varargin)

%MODIFY modify instance of STREAMobj to meet user-defined criteria
%
% Syntax
%
%     S2 = modify(S,pn,pv)
%     [S2,nalix] = ...
%
% Description
%
%     The function modify changes the geometry of a stream network to meet
%     different user-defined criteria. See
%
%     demo_modifystreamnet
%
%     for an (incomplete) overview of the function's scope. See input
%     arguments below for an overview on the functionality of the function.
%
% Input arguments
%
%     S     instance of stream object
%
%     Parameter name/value pairs (one pn/pv pair is required to run the
%     function)
%
%     'streamorder'   scalar or string
%     scalar: if you supply a scalar integer x as parameter value, S2
%     contains only streams with order x.
%     string: a string allows you to define a simple relational statement,
%     e.g. '>=3' returns the stream network with streams of streamorder
%     greater or equal than 3. Other commands may be '==5', '<5', etc.
%
%     'distance' scalar or [2x1] vector
%     lets you define a minimum and maximum distance from which the stream 
%     network starts and ends based on the distances given in S. A scalar is
%     interpreted as minimum upstream distance and a two element vector is
%     interpreted as minimum and maximum upstream distance to which the
%     stream network is cut.
%
%     'maxdsdistance' scalar
%     modifies the stream network that only the portion of the network is
%     retained that is within downstream distance of the the channelheads.
%
%     'upstreamto' logical GRIDobj, nal or linear index in GRIDobj
%     returns the stream network upstream to true pixels in the logical
%     raster of GRIDobj. Note that, if the grid contains linear features
%     (e.g. a fault), the line should be 4 connected. Use
%     bwmorph(I.Z,'diag') to establish 4-connectivity in a logical raster
%     I.
%
%     'downstreamto' logical GRIDobj, nal or linear index in GRIDobj
%     returns the stream network downstream to true pixels in the logical
%     raster of GRIDobj. Note that, if the grid contains linear features
%     (e.g. a fault), the foreground (true pixels) should be 4 connected.
%     Use bwmorph(I.Z,'diag') to establish 4-connectivity in a logical
%     raster I.
%
%     'rmconncomps' scalar
%     removes connected components (individual stream 'trees') of the
%     entire network with a maximum distance in map units less than the 
%     specified value.
%
%     'rmconncomps_ch' scalar
%     removes connected components (individual stream 'trees') that have
%     less or equal the number of channel heads 
%
%     'rmupstreamtoch' STREAMobj
%     removes nodes from S that are upstream to the channelheads of the
%     supplied STREAMobj.
%
%     'rmnodes' STREAMobj 
%     removes nodes in S that belong to another stream network S2.
%
%     'tributaryto' instance of STREAMobj
%     returns the stream network that is tributary to a stream (network) 
%
%     'tributaryto2' instance of STREAMobj
%     same as 'tributaryto' but tributaries include the pixel of the
%     receiving stream.
%
%     'lefttrib' instance of STREAMobj
%     returns the stream network that is tributary to a stream from the
%     left. For example, modify(S,'lefttrib',trunk(S)) selects the streams
%     in S that are left tributaries to the trunk stream of S.
%
%     'righttrib' instance of STREAMobj
%     returns the stream network that is tributary to a stream from the
%     right.
%
%     'shrinkfromtop' scalar distance
%     removes parts of the stream network that are within the specified
%     distance from the channelheads.
%
%     'clip' n*2 matrix or logical GRIDobj
%     retains those parts of the network that are inside the polygon
%     with vertices defined by the n*2 matrix with x values in the first
%     column and y values in the second column. The function automatically
%     closes the polygon if the first and the last row in the matrix
%     differ.
%
%     'interactive'  string
%        'polyselect': plots the stream network and starts a polygon tool 
%                      to select the stream network of interest.
%        'outletselect': select a basin based on a manually selected
%                      outlet.
%        'reachselect': select a reach based on two locations on the network
%        'channelheadselect': select a number of channel heads and derive 
%                      stream network from them.
%        'rectselect': like 'polyselect', but with a rectangle
%        'ellipseselect': like 'polyselect', but with an ellipse
%
%
% Output arguments
%
%     S2      modified stream network (class: STREAMobj)
%     nalix   index into node attribute list nal of S, so that nal2 =
%             nal(nalix)
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S  = STREAMobj(FD,flowacc(FD)>1000);
%
%     % Extract the stream network that is tributary to the main trunk
%     St = modify(S,'tributaryto',trunk(S));
%     plot(S)
%     hold on
%     plot(St,'r')
%     holf off
%
%     % Split the stream network into the network downslope and upslope to
%     % the 1200 m contour
%     C = dilate(DEM-1200,ones(3)) > 0 & erode(DEM-1200,ones(3)) < 0;
%     C.Z = bwmorph(C.Z,'skel',inf);
%     C.Z = bwmorph(C.Z,'diag');
%     Su  = modify(S,'upstreamto',C);
%     Sl  = modify(S,'downstreamto',C);
%     plot(Su)
%     hold on
%     plot(Sl,'r')
%     hold off
%
%     
% See also: STREAMobj, STREAMobj/trunk, STREAMobj/subgraph, 
%           demo_modifystreamnet
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. June, 2019

narginchk(3,3)

% Parse Inputs
p = inputParser;         
p.FunctionName = 'modify';
addRequired(p,'S',@(x) isa(x,'STREAMobj'));

addParamValue(p,'streamorder',[]);
addParamValue(p,'distance',[],@(x) isnumeric(x) && numel(x)<=2);
addParamValue(p,'maxdsdistance',[],@(x) isscalar(x) && x>0);
addParamValue(p,'interactive',[],@(x) ischar(validatestring(x,{'polyselect','reachselect',...
                                                               'channelheadselect','rectselect',...
                                                               'ellipseselect','outletselect'})));
addParamValue(p,'tributaryto',[],@(x) isa(x,'STREAMobj'));
addParamValue(p,'tributaryto2',[],@(x) isa(x,'STREAMobj'));
addParamValue(p,'righttrib',[],@(x) isa(x,'STREAMobj'));
addParamValue(p,'lefttrib',[],@(x) isa(x,'STREAMobj'));
addParamValue(p,'shrinkfromtop',[],@(x) isnumeric(x) && isscalar(x) && x>0);
addParamValue(p,'upstreamto',[],@(x) isa(x,'GRIDobj') || isnumeric(x)  || isnal(S,x));
addParamValue(p,'downstreamto',[],@(x) isa(x,'GRIDobj') || isnumeric(x) || isnal(S,x));
addParamValue(p,'rmconncomps',[],@(x) isnumeric(x) && x>0 && isscalar(x));
addParamValue(p,'rmconncomps_ch',[],@(x) isnumeric(x) && x>=0 && isscalar(x));
addParamValue(p,'rmupstreamtoch',[],@(x) isa(x,'STREAMobj'));
%addParamValue(p,'between',[],@(x) isnumeric(x) && size(x,2) == 2);
addParamValue(p,'fromch2IX',[]);
addParamValue(p,'rmnodes',[]);
addParamValue(p,'clip',[],@(x) (isnumeric(x) && size(x,2)==2 && size(x,1)>2) || isa(x,'GRIDobj'));
addParamValue(p,'nal',[],@(x) isnal(S,x));

parse(p,S,varargin{:});
S   = p.Results.S;


if ~isempty(p.Results.streamorder)
%% streamorder    
    so = streamorder(S);
    if isnumeric(p.Results.streamorder)
        validateattributes(p.Results.streamorder,{'numeric'},...
            {'scalar','integer','>',0,'<=',max(so)});
        I = so == p.Results.streamorder;
    else
        try
            % expecting relational operator
            sothres = nan;
            counter = 1;
            while isnan(sothres)
                sothres = str2double(p.Results.streamorder(counter+1:end));
                relop   = p.Results.streamorder(1:counter);
                counter = counter+1;
            end
        catch ME
            rethrow(ME)
        end

        relop = str2func(relop);
        I = relop(so,sothres);
        
    end


elseif ~isempty(p.Results.distance)
%% distance        
    d = S.distance;
    maxdist  = max(d);
    distrange = p.Results.distance;
    
    
    validateattributes(distrange,{'numeric'},...
        {'>=',0,'<=',maxdist});
    
    if isscalar(distrange)
        distrange(2) = maxdist;
    end
    
    if distrange(2) <= distrange(1)
        error('TopoToolbox:wronginput',...
            'The minimum distance must not be equal to or greater than the maximum distance');
    end
        
    
    I = d>=distrange(1) & d<=distrange(2);               % + norm([S.cellsize S.cellsize]);
    
elseif ~isempty(p.Results.maxdsdistance)
%% maximmum downstream distance    
    d = distance(S,'min_from_ch');
    I = d <= p.Results.maxdsdistance;
    
elseif ~isempty(p.Results.upstreamto)
%% upstream to    
    
    II = p.Results.upstreamto;
    
    if isa(II,'GRIDobj')
        
        II = getnal(S,II);
    elseif isnal(S,II)
        II = logical(II);        
    else   
        % II contains linear indices into the DEM from which S was derived
        [II,locb] = ismember(S.IXgrid,II);
%         if ~all(locb)
%             warning('TopoToolbox:modify','Some linear indices are not located on the stream network')            
%         end
    end
    
%     I = false(size(S.x));
    I = II;
    for r = numel(S.ix):-1:1
        I(S.ix(r)) = II(S.ixc(r)) || I(S.ixc(r))  || I(S.ix(r));
    end
% elseif ~isempty(p.Results.between)
%     %% between
%     M = sparse(S.ix,S.ixc,1,numel(S.x),numel(S.x));
%     G = graph(M,T);
%     ix = p.Results.between;
%     [II,ixs] = ismember(ix,S.IXgrid);
%     if any(II==0)
%         error('All cells must be on the network.')
%     end
%     % sort indices
%     [~,temp] = sort(S.distance(ixs),1,'ascend');
%     ixs = ixs()
% 
%     I  = logical(getnal(S));
%     for r = 1:size(ix,1)
%         [II,ixs] = ismember(S.IXgrid,ix(r,:));
%         
%         I(shortestpath(G,ixs))
%     end
    
elseif ~isempty(p.Results.rmupstreamtoch)
    
    S2 = p.Results.rmupstreamtoch;
    ch = streampoi(S2,'channelhead','logical');
    ch = nal2nal(S,S2,ch,false);
    
    I  = ch;
    
    for r = numel(S.ix):-1:1
        I(S.ix(r)) = I(S.ix(r)) || I(S.ixc(r));
    end
    
    I  = ~I;
    % retain channelheads in the list
    I(ch) = true;
    
        
elseif ~isempty(p.Results.downstreamto)
%% downstream to    
    
    II = p.Results.downstreamto;
    
    if isa(II,'GRIDobj')
        
        II = getnal(S,II);
    elseif isnal(S,II)
        II = logical(II);        
    else   
        % II contains linear indices into the DEM from which S was derived
        [II,locb] = ismember(S.IXgrid,II);
%         if ~all(locb)
%             warning('TopoToolbox:modify','Some linear indices are not located on the stream network')            
%         end
    end
    
    
    I = II;
    for r = 1:numel(S.ix)
        I(S.ixc(r)) = II(S.ix(r)) || I(S.ix(r)) || I(S.ixc(r));
    end

elseif ~isempty(p.Results.tributaryto) || ~isempty(p.Results.tributaryto2)
%% tributary to
    if ~isempty(p.Results.tributaryto)
        Strunk = p.Results.tributaryto;
    else
        Strunk = p.Results.tributaryto2;
    end
    
    II = ismember(S.IXgrid,Strunk.IXgrid);
    I  = false(size(S.x));
    
    for r = numel(S.ix):-1:1
        I(S.ix(r)) = (II(S.ixc(r)) || I(S.ixc(r))) && ~II(S.ix(r));
    end

    if ~isempty(p.Results.tributaryto2)
        % add trunk stream pixels
        I(S.ixc(II(S.ixc) & ~II(S.ix))) = true;
    end

elseif ~isempty(p.Results.righttrib) || ~isempty(p.Results.lefttrib)
%% select tributaries from a specified direction    
    % calculate directions
    direc = tribdir(S);
    
    if ~isempty(p.Results.righttrib)
        St = p.Results.righttrib;
        val = 1;
    else
        St = p.Results.lefttrib;
        val = -1;
    end
    II = STREAMobj2GRIDobj(St);
    ii = II.Z(S.IXgrid(S.ixc)) & ~II.Z(S.IXgrid(S.ix)) & (direc(S.ix) == val);
    
    I = false(size(S.x));
    for r = numel(S.ix):-1:1
        I(S.ix(r)) = I(S.ixc(r)) || ii(r);
    end
    
elseif ~isempty(p.Results.shrinkfromtop)
%% shrink from top
    d = distance(S,'max_from_ch');
    I = d > p.Results.shrinkfromtop;
elseif ~isempty(p.Results.fromch2IX)
    ch = streampoi(S,'channelhead','logical');
    en = ismember(S.IXgrid,p.Results.fromch2IX);
    
    pp = en;
    for r = numel(S.ix):-1:1
        if pp(S.ixc(r))
            pp(S.ix(r)) = true;
        end
    end
     
    for r = 1:numel(S.ix)
        if ch(S.ix(r)) && ~en(S.ixc(r)) && pp(S.ix(r))
            ch(S.ixc(r)) = true;
        end
    end
    
    I = ch;
    
elseif ~isempty(p.Results.rmconncomps)
%% remove connected conn comps
    cc = conncomps(S);
    d  = S.distance;
    md = find(accumarray(cc,d,[max(cc) 1],@max) > p.Results.rmconncomps);
    I = ismember(cc,md);

elseif ~isempty(p.Results.rmconncomps_ch)
%% remove connected conn comps
    cc = conncomps(S);
    d  = distance(S,'nr_of_ch');
    md = find(accumarray(cc,d,[max(cc) 1],@max) > p.Results.rmconncomps_ch);
    I = ismember(cc,md);
    
elseif ~isempty(p.Results.rmnodes)
    % check if empty S is supplied
    if isempty(p.Results.rmnodes.IXgrid)
        I = getnal(S) == 0;
    else
        I = ~ismember(S.IXgrid,p.Results.rmnodes.IXgrid);
    end

elseif ~isempty(p.Results.nal)
    I = p.Results.nal;

elseif ~isempty(p.Results.clip)
    if isa(p.Results.clip,'GRIDobj')
        mask = p.Results.clip;
        I    = mask.Z(S.IXgrid) > 0;
    else
        pos = p.Results.clip;
        if ~isequal(pos(end,:),pos(1,:))
            pos(end+1,:) = pos(1,:);
        end
        I = inpolygon(S.x,S.y,pos(:,1),pos(:,2));
    end

elseif ~isempty(p.Results.interactive)
%% interactive    
%     figure
    plot(S,'k'); axis equal
    ax = gca;
    % expand axes
    lims = axis;
    xlim(ax,[lims(1)-(lims(2)-lims(1))/20 lims(2)+(lims(2)-lims(1))/20])
    ylim(ax,[lims(3)-(lims(4)-lims(3))/20 lims(4)+(lims(4)-lims(3))/20])
    
    meth = validatestring(p.Results.interactive,...
        {'polyselect','reachselect','channelheadselect',...
         'rectselect','ellipseselect','outletselect'});
    switch meth
        case 'outletselect'
            title('map outlet and press key to finalize')
            hold on
            htemp = plot(ax,[],[]);
            hp = impoint('PositionConstraintFcn',@getnearest);
            addNewPositionCallback(hp,@drawbasin);
            setPosition(hp,getPosition(hp));
            set(gcf,'WindowKeyPressFcn',@(k,l) uiresume);
            pos = wait(hp);
%             pos = getPosition(hp);
%             
            hold off
            delete(hp);
            I = S.x==pos(1) & S.y==pos(2);
            for r = numel(S.ixc):-1:1
                I(S.ix(r)) = I(S.ixc(r)) || I(S.ix(r));
            end
            
        
        case 'channelheadselect'
            title('map channel heads and enter any key to finalize')
            hold on
            xy = streampoi(S,'channelhead','xy');
            scatter(xy(:,1),xy(:,2));
            hold off
            set(gcf,'WindowKeyPressFcn',@(k,l) uiresume);
            xys   = [];
            while true
                try
                    hpstart = impoint('PositionConstraintFcn',@getnearestchanhead);                    
                    setColor(hpstart,[0 1 0])
                    setPosition(hpstart,getPosition(hpstart))
                    xys = [xys; getPosition(hpstart)]; %#ok<AGROW>
                catch
                    break
                end
            end
            
            if ~isempty(xys)
                I = ismember([S.x S.y],xys,'rows');
            else
                error('You must choose at least one channelhead')
            end
            
            for r = 1:numel(S.ix)
                I(S.ixc(r)) = I(S.ix(r)) || I(S.ixc(r));
            end
               
        case {'polyselect', 'rectselect', 'ellipseselect'}
            switch meth
                case 'polyselect'
                    title('create a polygon and double-click to finalize')
                    try
                    hp = drawpolygon;
                    pos = hp.Position;
                    catch
                    hp = impoly;
                    pos = wait(hp);
                    end
                    pos(end+1,:) = pos(1,:);
                case 'ellipseselect'
                    title('create a ellipse and double-click to finalize')
                    try
                    hp = drawellipse;
                    pos = hp.Vertices;
                    catch
                    hp = imellipse;
                    pos = wait(hp);
                    pos = getVertices(hp);
                    end
                    
                    pos(end+1,:) = pos(1,:);
                case 'rectselect'
                    title('create a rectangle and double-click to finalize')
                    hp = imrect;
                    pos = wait(hp);
                    pos = [pos(1)         pos(2); ...
                           pos(1)+pos(3)  pos(2); ...
                           pos(1)+pos(3)  pos(2)+pos(4);...
                           pos(1)         pos(2)+pos(4);...
                           pos(1)         pos(2)];
            end

            I = inpolygon(S.x,S.y,pos(:,1),pos(:,2));
            hold on
            plot(pos(:,1),pos(:,2),'-');
            hold off
            delete(hp);
        case 'reachselect'
            ix = [];
            hpath = [];
            ixpath = [];
            ixend  = 0;
            ixcix = zeros(numel(S.x),1);
            ixcix(S.ix) = 1:numel(S.ix);
            title('set upper (green) reach location')
            hpstart = impoint('PositionConstraintFcn',@getnearest);
            setColor(hpstart,[0 1 0])
            idstart = addNewPositionCallback(hpstart,@drawpath);
            setPosition(hpstart,getPosition(hpstart))
            
            title('set lower (red) reach location')
            hpend   = impoint('PositionConstraintFcn',@getnearestend);            
            setColor(hpend,  [1 0 0])
            idend = addNewPositionCallback(hpend,@drawpath);
            setPosition(hpend,getPosition(hpend))
            
            title('move reach locations and press any key to extract reach')
            set(gcf,'WindowKeyPressFcn',@(k,l) uiresume);
            uiwait
            
            I = false(size(S.x));
            I(ixpath) = true;
            
            delete(hpstart)
            delete(hpend);


    end
end

%% clean up
if exist('I','var')

    Sout = subgraph(S,I);
    Sout = clean(Sout);
    
    if nargout == 2
        [~,nalix] = ismember(Sout.IXgrid,S.IXgrid); 
    end

end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function posn = getnearest(pos)
    % GETNEAREST  snap to nearest location on stream network
    [~,ix] = min((S.x-pos(1)).^2 + (S.y-pos(2)).^2);
    posn= [S.x(ix) S.y(ix)];

end

function posn = getnearestchanhead(pos)
    % GETNEAREST  snap to nearest location on stream network
    [~,ix] = min((xy(:,1)-pos(1)).^2 + (xy(:,2)-pos(2)).^2);
    posn= xy(ix,:);

end

function posn = getnearestend(pos)
    % GETNEAREST  snap to nearest location on stream network
        [~,ixend] = min((S.x-pos(1)).^2 + (S.y-pos(2)).^2);
        posn= [S.x(ixend) S.y(ixend)];

end

function drawpath(pos)
    % DRAWPATH   
    
    % bring hpstart and hpend to top

    ixcc = ix;
    ixpath = ix;
    while ixcix(ixcc) ~= 0 && ixcc ~= ixend %S.ixc(ixcix(ixcc)) ~= ixend;
        ixpath(end+1) = S.ixc(ixcix(ixcc));
        ixcc = ixpath(end);
    end

    if ishandle(hpath)
        set(hpath,'Xdata',S.x(ixpath),'Ydata',S.y(ixpath));
    else
        hold on
        hpath = plot(S.x(ixpath),S.y(ixpath),'r','LineWidth',1.5);
        hold off
    end
end 

    function II = drawbasin(pos)
        II = S.x==pos(1) & S.y==pos(2);
        IX = S.IXgrid(II);
        Stemp = modify(S,'upstreamto',IX);
        delete(htemp)
        try
            htemp = plot(Stemp,'r');
        end
        
    end
        
end
                        
        
    