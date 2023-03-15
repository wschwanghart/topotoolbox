function OUT = flowdistance(FD,varargin) 

%FLOWDISTANCE flow distance in upstream and downstream direction
%
% Syntax
%
%     D = flowdistance(FD)
%     D = flowdistance(FD,GRIDobj)
%     D = flowdistance(FD,S)
%     D = flowdistance(FD,IX)
%     D = flowdistance(FD,x,y)
%     D = flowdistance(...,direction)
%
% Description
%
%     flowdistance calculates the horizontal distance along the flow
%     network in either upstream (default) or downstream direction.
%     Downstream flow distance is the maximum distance along the drainage
%     network (flow length). If only a flow direction object (FLOWobj) is
%     supplied to the function, flowdistance calculates the distance from
%     outlets and ridges in upstream and downstream direction,
%     respectively. If seed locations for the distance transform are
%     supplied (e.g. logical GRIDobj, STREAMobj, linear index or coordinate 
%     pairs) the distance is calculated from the seed locations. When using
%     a STREAMobj as input, usually upstream distance makes sense only.
%
% Input arguments
%
%     FD          flow direction (FLOWobj)
%     GRIDobj     logical grid (GRIDobj)
%     S           instance of STREAMobj
%     IX          linear index of seeds
%     x,y         x- and y-coordinate vectors
%     direction   'upstream' (default) or 'downstream' (maximum flow path
%                 length). By default, 'downstream' is 'maxdownstream'
%                 which means that the maximum downstream distance is
%                 calculated in downstream direction. 'mindownstream'
%                 calculates the minimum downstream distance.
%
% Output arguments
%
%     D           distance grid (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     D = flowdistance(FD);
%     imageschs(DEM,D)
%
% See also: FLOWobj, FLOWobj/flowacc, GRIDobj
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 30. November, 2021


% 4/3/2016: the function now makes copies of FD.ix and FD.ixc (see 
% FLOWobj/flowacc
% 11/30/2021: update of help and making the mindownstream option visible.

%% check input arguments

narginchk(1,4)

% check if last input argument is a string 
if ~isempty(varargin) && ischar(varargin{end})
    direction = validatestring(varargin{end},{'upstream','downstream','maxdownstream','mindownstream'});
    nrargsin = nargin - 1;
else
    direction = 'upstream';
    nrargsin = nargin;
end



if nrargsin == 2
    % SEED pixels are either supplied as logical matrix, GRIDobj, or linear
    % index, or STREAMobj
    SEED = varargin{1};
    isGRIDobj = isa(SEED,'GRIDobj');   
    if (islogical(SEED) || isGRIDobj)
        validatealignment(FD,SEED);
        if isGRIDobj
            SEED = SEED.Z;
        end
    elseif isa(SEED,'STREAMobj')
        ix = SEED.IXgrid;
        SEED = false(FD.size);
        SEED(ix) = true;
        
        
    else
        % SEED is supposed to be supplied as linear index in the GRIDobj
        ix   = varargin{1};
        SEED = false(FD.size);
        SEED(ix) = true;
    end
elseif nrargsin == 3
    % SEED pixels are supplied as coordinate pairs
    ix   = coord2ind(FD,varargin{1},varargin{2});
    SEED = false(FD.size);
    SEED(ix) = true;
end

if ~strcmpi(FD.type,'single')
    error('TopoToolbox:FLOWobj','flowdistance requires flow distance type to be single')
end

%% Do calculation
ixtemp  = FD.ix;
ixctemp = FD.ixc;

% Is FD geographic? 
% if ~isgeographic(FD)
    DIST = getdistance(ixtemp,ixctemp,FD.size,FD.cellsize,'single');
% else
%     [LON,LAT] = getcoordinates(GRIDobj(FD,'logical'),'matrix');
%     DIST = sph_distance(LAT(ixtemp),LON(ixtemp),LAT(ixctemp),LON(ixctemp),FD.georef.gcs);
%     clear LON LAT
% end

cl   = class(DIST);
switch direction
    case 'upstream'
        %% Upstream distance calculation
        if nrargsin == 1
            D     = zeros(FD.size,cl);
            start = numel(ixtemp);
            for r = start:-1:1
                D(ixtemp(r)) = D(ixctemp(r))+DIST(r);
            end
        else
            D     = inf(FD.size,cl);
            D(SEED) = 0;
            start = find(SEED(ixctemp),1,'last');
            for r = start:-1:1
                D(ixtemp(r)) = min(D(ixctemp(r))+DIST(r),D(ixtemp(r)));
            end
            
            D(isinf(D)) = nan;
        end
    case {'downstream' 'maxdownstream'}
        
        %% Downstream distance calculation
        
        if nrargsin == 1
            D = zeros(FD.size,cl);
        else
            D = -inf(FD.size,cl);
            D(SEED) = 0;
        end
        
        for r = 1:numel(FD.ix)
            D(FD.ixc(r)) = max(D(FD.ix(r))+DIST(r),D(FD.ixc(r)));
        end
        
        if nrargsin >= 2
            D(isinf(D)) = nan;
        end
    case 'mindownstream'
              
        
        D = inf(FD.size,cl);
        D(SEED) = 0;
        
        for r = 1:numel(ixtemp)
            D(ixctemp(r)) = min(D(ixtemp(r))+DIST(r),D(ixctemp(r)));
        end
end


%% Prepare Output
OUT = GRIDobj(FD);
% write output to GRIDobj
OUT.Z = D;
OUT.zunit = '';
OUT.name  = [direction ' flow distance'];


