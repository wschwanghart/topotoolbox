function [OUT,varargout] = drainagebasins(FD,varargin)

%DRAINAGEBASINS drainage basin delineation/catchments
%
% Syntax
%
%     L = drainagebasins(FD)
%     L = drainagebasins(FD,IX)
%     L = drainagebasins(FD,S)
%     L = drainagebasins(FD,x,y)
%     L = drainagebasins(FD,SO,order)
%     [L,outlets] = ...
%     [L,x,y] = ...
%
% Description
%
%     drainagebasins determines drainage basin affiliation of cells in 
%     a digital elevation model (DEM) based on flow direction (FLOWobj).
%     When drainagebasins is supplied with the flow direction matrix, all
%     cells have a clear affiliation to a drainage basin.
%
%     drainagebasins(FD) calculates the drainage basins from all drainage
%     network outlets.
%
%     drainagebasins(FD,IX) and drainagebasins(FD,x,y) lets you specify
%     outlets indicated by the linear index vector (IX) or coordinate
%     vectors (x,y) to derive drainage basins.
%
%     drainagebasins(FD,S) takes the outlets of the stream network in the
%     STREAMobj S.
%
%     drainagebasins(FD,SO,order) takes a streamorder grid (SO) calculated
%     by the function FLOWobj/streamorder and a scalar that indicates the
%     stream order for which drainage basins are to be derived. order = 1 
%     for instance, outputs all first order drainage basins. Note that
%     drainage basins are not calculated for stream links of specified 
%     order that are attached to the grids edges. 
%
% Input
%
%     FD        flow direction object (FlowDirObj)
%     IX        linear index indicating position of user defined drainage 
%               basin outlets
%     S         STREAMobj. Takes the outlets of the stream network S.
%     x,y       coordinate pairs of user defined drainage basin outlets
%     SO        stream order raster as produced by the function
%               FLOWobj/streamorder (GRIDobj)
%     order     stream order (scalar integer)
% 
% Output
%
%     L         drainage basin grid (GRIDobj, uint32)
%     outlets   linear indices of drainage basin outlet cells
%     x,y       coordinate pairs of drainage basin outlet cells
%
% Example 1
%
%     % get all drainagebasins
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     D   = drainagebasins(FD);
%     imageschs(DEM,D)
%
% Example 2
%
%     % get all drainage basins with streamorder 3
%     S = streamorder(FD,flowacc(FD)>100);
%     D = drainagebasins(FD,S,3);
%     imageschs(DEM,D)
%
% Example 3: map values to drainage basins (e.g. basin wide erosion rates).
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S  = STREAMobj(FD,'minarea',1000);
%     IX = randlocs(S,10);
%     erosrate = rand(size(IX));
%     D  = drainagebasins(FD,IX);
%     E = GRIDobj(DEM)*nan;
%     E.Z(D.Z~=0) = erosrate(D.Z(D.Z~=0));
%     imageschs(DEM,E)
%
%
% See also: FLOWobj, GRIDobj/snap2stream, FLOWobj/streamorder,
%           GRIDobj/shufflelabel
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 3. December, 2018



% 4/3/2016: the function now makes copies of FD.ix and FD.ixc (see 
% FLOWobj/flowacc

% 30/9/2016: bug removed when called with 2 or 3 outputs

% 3/12/2018: new example and some additional error checking

narginchk(1,3);

switch lower(FD.type)
    case {'multi','dinf'}
    error('TopoToolbox:drainagebasins',...
        'Drainage basins are not defined for divergent flows');
end

% create temporary indices to gain speed with R2015b and later
ixtemp  = FD.ix;
ixctemp = FD.ixc;


% If the function is called with only one output and input argument, and if
% the mex-function exists, then the mex-function will be called
if exist(['drainagebasins_mex.' mexext],'file')==3 && ...
   nargout==1 && nargin == 1
    % Use mex-file
    D = drainagebasins_mex(ixtemp,ixctemp,FD.size);

elseif nargin == 1    
    % Don't use mex-file
    DBcounter = 0;
    D = zeros(FD.size,'uint32');

    for r = numel(ixtemp):-1:1
        if D(ixctemp(r)) == 0
            DBcounter = DBcounter+1;
            D(ixctemp(r)) = DBcounter;
            outlets(DBcounter) = ixctemp(r);
        end
        D(ixtemp(r)) = D(ixctemp(r));
    end
    
    outlets = double(outlets(:));
    
elseif nargin > 1
    % ix,x and y, or Stream order grid and stream order are supplied
    if nargin == 2
        if isa(varargin{1},'STREAMobj')
            IX = streampoi(varargin{1},'outlets','ix');
        else
            % IX is supplied
            IX = varargin{1}; 
            validateattributes(IX,{'numeric'},{'vector','integer'})
        end
    else
        if isa(varargin{1},'GRIDobj')
            S = varargin{1};
            % find seed pixels for specific stream order
            S.Z(S.Z ~= varargin{2}) = 0;            
            I = (S.Z(ixtemp) & ~S.Z(ixctemp));
            IX = ixtemp(I);
        else
            % SEED pixels are supplied as coordinate pairs
            IX   = coord2ind(FD,varargin{1},varargin{2});
        end
    end
    
    if any(IX<1) || any(IX>prod(FD.size))
        error('TopoToolbox:drainagebasins',...
            'Some of the locations are outside the grid boundaries.')
    end
 
    D = zeros(FD.size,'uint32');
    D(IX) = cast(1:numel(IX),'uint32');
        
    for r = numel(ixtemp):-1:1
        if D(ixctemp(r)) ~= 0 && D(ixtemp(r))==0
            D(ixtemp(r)) = D(ixctemp(r));
        end
    end
    
    outlets = IX;
end

%% Prepare Output
OUT = copy2GRIDobj(FD);
% write output to GRIDobj
OUT.Z = D;
OUT.zunit = '';
OUT.name  = 'drainage basins';

if nargout == 2
    varargout{1} = outlets;
elseif nargout == 3
    [x,y] = ind2coord(FD,outlets);
    varargout{1} = x;
    varargout{2} = y;
end


    