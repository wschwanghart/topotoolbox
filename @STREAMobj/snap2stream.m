function [xn,yn,IX,D,MP] = snap2stream(S,x,y,varargin)

%SNAP2STREAM snap locations to nearest stream location
%
% Syntax
%
%     [xn,yn] = snap2stream(S,x,y)
%     [xn,yn] = snap2stream(S,x,y,pn,pv,...)
%     [xn,yn] = snap2stream(S,lat,lon,'inputislatlon',true,pn,pv,...)
%     [xn,yn,IX,res,MP] = snap2stream(...)
%     ... = snap2stream(...,'nalarea',an,'pointarea',ap,'maxdist',d)
%
% Description
%
%     Often gauge locations are not precisely associated with 
%     streams derived from a digital elevation model. Calculation of 
%     drainage basins may thus return wrong basin delineations. One method 
%     is to manually adjust their location to a derived stream raster.
%     snap2stream does this automatically by finding identifying the
%     shortest distance to a stream. 
%
%     Sometimes, however, assigning nearest locations may lead to snapping
%     to the wrong locations (e.g. tributaries instead to the trunk river).
%     Providing additional argument pairs 'nalarea' and 'pointarea' offer
%     additional means to constrain the snapping. The unit of the areas
%     should be in the same units as the horizontal coordinates (e.g. m). 
%
% Input
% 
%     S      instance of STREAMobj
%     x,y    coordinate vectors
%     
%     Parameter name/value pairs {default}
%
%     'snapto'   string 
%     constrain snapping to network features. Features can be {'all'},
%     'confluence', 'outlet' or 'channelhead'
%
%     'streamorder'   scalar or string
%     scalar: if you supply a scalar integer s as parameter value,
%     locations will be snapped only to streams with order s. string: a
%     string allows you to define a simple relational statement, e.g. '>=3'
%     will snap only to streams of streamorder greater or equal than 3.
%     Other commands may be '==5', '<5', etc.
%
%     'maxdist'  scalar
%     maximum snapping distance {inf}. 
%
%     'inputislatlon'   true or {false}
%     determines whether input coordinates are supplied in a geographic
%     coordinate system (wgs84). If set to true, snap2stream projects the
%     geographic coordinates to the same coordinate system as S. This
%     requires that S has a projected coordinate system and the mapping
%     toolbox. 
%
%     'alongflow' [] or FLOWobj
%     if a FLOWobj FD is provided, then points are moved along the flow
%     network FD to the snap location. (nalarea and pointarea cannot be
%     used together with this option)
%
%     'nalarea'   []
%     Node-attribute list with upslope area values for each node in the
%     river network (e.g. getnal(S,flowacc(FD)*DEM.cellsize^2))
%
%     'pointarea'  []
%     Upslope area for each location in the same unit as in 'nalarea' 
%
%     'plot'  true or {false}
%     plot results
%
%
% Output arguments
%
%     xn,yn    coordinate vectors of snapped locations
%     IX       linear index into grid from which S was derived
%     res      residual (eucl. distance between locations)
%     MP       map struct that can be exported to a shapefile
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,A>1000);
%     IX  = randperm(prod(DEM.size),20);
%     [x,y] = ind2coord(DEM,IX);
%     [xn,yn,IX,res,d,MP] = snap2stream(S,x,y,'alongflow',FD,'plot',true);
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 12. January, 2021


% Parse Inputs
p = inputParser;         
p.FunctionName = 'snap2stream';

snaptomethods = {'all','confluence','outlet','channelhead'};

addParamValue(p,'snapto','all',@(x) ischar(validatestring(x,snaptomethods)));
addParamValue(p,'maxdist',inf,@(x) isscalar(x) && x>0);
addParamValue(p,'streamorder',[]);
addParamValue(p,'inputislatlon',false,@(x) isscalar(x));
addParamValue(p,'plot',false,@(x) isscalar(x));
addParamValue(p,'nalarea',[],@(x) isnal(S,x));
addParamValue(p,'pointarea',[],@(x) numel(x) == numel(x));
addParamValue(p,'alongflow',[],@(x) isa(x,'FLOWobj') || isempty(x));

parse(p,varargin{:});

if ~isempty(p.Results.nalarea)
    if isempty(p.Results.pointarea) || isempty(p.Results.maxdist)
        error('TopoToolbox:snap2stream',...
            ['snap2stream requires the arguments pointarea and maxdist\newline'...
            'if nalarea is provided']);
    end
end
        

snapto = validatestring(p.Results.snapto,{'all','confluence','outlet','channelhead'});
I      = true(size(S.x));

if p.Results.inputislatlon
    if isempty(S.georef)
        error('TopoToolbox:georeferencing',...
            ['Projection of stream network unknown. Cannot convert ' ...
             'geographic coordinates.']);
    end
    [x,y] = mfwdtran(S.georef.mstruct,x,y);
end

% evaluate optional arguments to restrict search to specific stream
% locations only

switch snapto
    case 'all'
        % do nothing
    otherwise
        nrc = numel(S.x);
        M = sparse(S.ix,S.ixc,true,nrc,nrc);
        
        switch snapto
            case 'confluence'
                I = I & (sum(M,1)>1)';
            case 'outlet'
                I = I & (sum(M,2) == 0);
            case 'channelhead'
                I = I & (sum(M,1) == 0)';
        end
end
        
if ~isempty(p.Results.streamorder)
    
    so = streamorder(S);
    if isnumeric(p.Results.streamorder)
        validateattributes(p.Results.streamorder,{'numeric'},...
            {'scalar','integer','>',0,'<=',max(so)});
        I = (so == p.Results.streamorder) & I;
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
        I = relop(so,sothres) & I;
        
    end
end

% find nearest stream location

xn = S.x(I);
yn = S.y(I);
IX = S.IXgrid(I);

xn = xn(:);
yn = yn(:);
IX = IX(:);

x  = x(:);
y  = y(:);

if isempty(p.Results.nalarea)
    % No area given
    
    if isempty(p.Results.alongflow)
        % fast and easy using Statistics toolbox
        try
            [IDX,D] = knnsearch([xn yn],[x y]);
        catch %#ok<CTCH>
            % aargh, not available... dirty implementation
            D = nan(numel(x),1);
            IDX = nan(numel(x),1);
            
            for r = 1:numel(x)
                [D(r), IDX(r)] = min(sqrt((xn-x(r)).^2 + (yn-y(r)).^2));
            end
        end
        xn = xn(IDX);
        yn = yn(IDX);
        IX = IX(IDX);
    
    else
        xn = nan(size(x));
        yn = nan(size(y));
        
        FD = p.Results.alongflow;
        % Get drainage basin of each stream pixel
        [DB,IXOUTLET] = drainagebasins(FD,IX);
        IXO = coord2ind(DB,x,y);
        % those outside the grid borders are invalid and are not snapped
        invalid = isnan(IXO);
        % there may be points in areas with no drainage basin, 
        % which we set to invalid, too
        invalid(~invalid) = DB.Z(IXO(~invalid))==0;
        xn(invalid) = x(invalid);
        yn(invalid) = y(invalid);
        IX = nan(size(xn));
        
        IXOUTLET = IXOUTLET(DB.Z(IXO(~invalid)));
        IX(~invalid) = IXOUTLET;
        [xn(~invalid),yn(~invalid)] = ind2coord(DB,IX(~invalid));
        
        DIST = flowdistance(FD,IX(~invalid));
        D    = inf(size(xn));
        D(~invalid) = DIST.Z(IXO(~invalid));
    end
    
else
    
    ap = p.Results.pointarea(:);
    an = p.Results.nalarea(I);
    
    xnal = xn;
    ynal = yn;
    ixnal = IX;
    [ix,CD] = rangesearch([xnal,ynal],[x y],p.Results.maxdist);
    xn = nan(size(ap));
    yn = nan(size(ap));
    IX = nan(size(ap));
    D  = nan(size(ap));
    for r = 1:numel(CD)
        if isempty(ix{r})
            continue
        end
        areadiff = abs(ap(r) - an(ix{r}));
        cost = sqrt(areadiff(:)) + CD{r}(:);
        [~,ii] = min(cost);
        xn(r) = xnal(ix{r}(ii));
        yn(r) = ynal(ix{r}(ii));
        IX(r) = ixnal(ix{r}(ii));
        D(r)  = CD{r}(ii);
    end
    
end



% apply maximum distance
if ~isinf(p.Results.maxdist)
    ID = D <= p.Results.maxdist;
else
    ID = true(numel(x),1);
end

if nargout == 5
    % create output mapstruct
    MP = struct('Geometry',{'Point'},...
        'X',num2cell(xn),...
        'Y',num2cell(yn),...
        'IX',num2cell(IX),...
        'residual',num2cell(D),...
        'inrange',num2cell(double(ID)),...
        'Xold',num2cell(x),...
        'Yold',num2cell(y));
end



if ~isinf(p.Results.maxdist)
    ID = ~ID;
    xn(ID) = nan;
    yn(ID) = nan;
    IX(ID) = nan;
    ID = ~ID;
end

% plot results
if p.Results.plot
    plot(S)
    hold on
    plot(x,y,'sk')
    plot(xn(ID),yn(ID),'xr')
    plot([xn(ID) x(ID)]',[yn(ID) y(ID)]','k');
end
end




