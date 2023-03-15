function OUT = dependencemap(FD,varargin)

%DEPENDENCEMAP upslope area for specific locations in a DEM
%
% Syntax
%
%     I = dependencemap(FD,L)
%     I = dependencemap(FD,ix)
%     I = dependencemap(FD,x,y)
%
% Description
%
%     dependencemap returns a GRIDobj with true values masking the upslope
%     part of the digital elevation model that contributes to the specified
%     area in the region in L or pixels with the linear index ix or
%     coordinates x,y.
%
% Input
%
%     FD        flow direction object (FlowDirObj)
%     L         logical grid (GRIDobj)
%     ix        linear index into a grid
%     x,y       coordinate pairs
%
% Output
%
%     I         logical influence grid (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     I   = GRIDobj(DEM,'logical');
%     I.Z(300:500,300:500) = true;
%     FD = FLOWobj(DEM,'preprocess','c');
%     D  = dependencemap(FD,I);
%     imageschs(DEM,I+D)
%
% See also: FLOWobj, FLOWobj/INFLUENCEMAP
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017


% 4/3/2016: the function now makes copies of FD.ix and FD.ixc (see 
% FLOWobj/flowacc


%% check input arguments
narginchk(1,3)
if nargin == 2;
    % SEED pixels are either supplied as logical matrix, GRIDobj, or linear
    % index
    SEED = varargin{1};
    isGRIDobj = isa(SEED,'GRIDobj');   
    if (islogical(SEED) || isGRIDobj)
        validatealignment(FD,SEED);
        if isGRIDobj
            SEED = SEED.Z;
        end
    else
        % SEED is supposed to be supplied as linear index in the GRIDobj
        ix   = varargin{1};
        SEED = false(FD.size);
        ix   = round(ix);
        if any(ix <= 0 | ix > prod(FD.size))
            error('TopoToolbox:WrongInput',...
            ['Linear indices must not be less or equal to zero \n' ...
             'or larger than ' double2str(prod(FD.size)) '.']);
        end
        
        SEED(ix) = true;
    end
elseif nargin == 3
    % SEED pixels are supplied as coordinate pairs
    ix   = coord2ind(FD,varargin{1},varargin{2});
    SEED = false(FD.size);
    SEED(ix) = true;
end


% this is a very crude and slow implementation of a graph traversal
% algorithm since all nodes are visited

if ~(exist(['dependencemap_mex.' mexext],'file') == 3)
    % m implementation
    ixtemp  = FD.ix;
    ixctemp = FD.ixc;
    for r = numel(ixtemp):-1:1
        SEED(ixtemp(r)) = SEED(ixtemp(r)) || SEED(ixctemp(r));
    end
    
else
    % mex implementation
    SEED = dependencemap_mex(FD.ix,FD.ixc,FD.size,find(SEED));    
end


%% Prepare Output
% empty GRIDobj
OUT = copy2GRIDobj(FD);
% write output to GRIDobj
OUT.Z = SEED;
OUT.zunit = 'logical';
OUT.name  = 'dependence map';


end








