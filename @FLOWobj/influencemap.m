function OUT = influencemap(FD,varargin)

%INFLUENCEMAP downslope area for specific locations in a digital elevation model
%
% Syntax
%
%     I = influencemap(FD,L)
%     I = influencemap(FD,ix)
%     I = influencemap(FD,x,y)
%
% Description
%
%     influencemap returns a GRIDobj with true values for those pixels that
%     are downstream of the locations in L, ix, or x and y.
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
%     FD = FLOWobj(DEM,'preprocess','c');
%     I  = influencemap(FD,540261);
%     imageschs(DEM,dilate(I,ones(5)));
%
%
% See also: FLOWobj, FLOWobj/DEPENDENCEMAP
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. March, 2016


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
        if isGRIDobj;
            SEED = SEED.Z;
        end
    else
        % SEED is supposed to be supplied as linear index in the GRIDobj
        ix   = varargin{1};
        SEED = false(FD.size);
        SEED(ix) = true;
    end
elseif nargin == 3;
    % SEED pixels are supplied as coordinate pairs
    ix   = coord2ind(FD,varargin{1},varargin{2});
    SEED = false(FD.size);
    SEED(ix) = true;
end

%% Do calculation
ixtemp  = FD.ix;
ixctemp = FD.ixc;
for r = 1:numel(FD.ix);
    SEED(ixctemp(r)) = SEED(ixtemp(r)) || SEED(ixctemp(r));
end

%% Prepare Output
% empty GRIDobj
OUT = copy2GRIDobj(FD);
% write output to GRIDobj
OUT.Z = SEED;
OUT.zunit = 'logical';
OUT.name  = 'influence map';
