function OUT = localtopography(DEM,varargin)

%LOCALTOPOGRAPHY Local topography
%
% Syntax
%
%     H = localtopography(DEM)
%     H = localtopography(DEM,radius)
%     H = localtopography(DEM,radius,pn,pv,...)
%
% Description
%
%     localtopography quantifies local relief, e.g. the elevation range
%     within a specific radius. localtopography may take a while to
%     evaluate large DEMs with large kernels. You may speed up calculations
%     by adjusting the parameter 'N' for calculating the max, min or range
%     filter. 
%
%     localtopography uses symmetric boundary padding. NaNs
%     are inpainted by nearest neighbor interpolation prior to the
%     calculation (affects only 'mean', 'median', 'prctile', 'std'). For
%     'max', 'min' and 'range', localtopography adopts the padding behavior
%     of the function imdilate and imerode.
%
%
% Input arguments
%
%     DEM    digital elevation model (GRIDobj)
%     radius radius of the moving window filter in map units. The default
%            value is 5000 (m).
%
% Parameter Name/Value (pn,pv) pairs
%
%     'type'   'range' (default), 'max', 'min', 'mean', 'median',
%              'prctile', 'std' (standard deviation)
%     'prc'    scalar between 0 and 100 [%]. Only applicable if 'prctile'
%              is chosen as 'type'.
%     'N'      enable speed improvement for 'max', 'min' and 'range'.
%              N must be 0, 4, 6,or 8 (default). When N is greater than 0,
%              the disk-shaped structuring element is approximated by a
%              sequence of N periodic-line structuring elements. When N
%              equals 0, no approximation is used, and the structuring
%              element members consist of all pixels whose centers are no
%              greater than R away from the origin. Choosing N entails a
%              trade-off between how well the shape of the disc-shaped
%              moving window is approximated and computational efficiency.
%              If you choose a large radius, then N=0 can take a while to 
%              evaluate.
%
% Output arguments
%
%     H      local topography grid (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     H = localtopography(DEM,500);
%     imageschs(DEM,H)
%
% See also: IMDILATE, IMERODE
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 28. January, 2013


narginchk(1,inf)

p = inputParser;
p.FunctionName = 'GRIDobj/localtopography';
expectedTypes = {'range','max','min','mean','median','prctile','std'};

addRequired(p,'DEM',@(x) issparse(x) || isa(x,'GRIDobj'));
addOptional(p,'radius',5000,@(x) isscalar(x) && x>DEM.cellsize);

addParamValue(p,'type','range',@(x) ischar(validatestring(x,expectedTypes)));
addParamValue(p,'N',8,@(x) isscalar(x) && ismember(x,[0 4 6 8]));
addParamValue(p,'thin',1,@(x) x>0.1 && x<=1);
addParamValue(p,'prc',90,@(x) x>0 && x<100);

parse(p,DEM,varargin{:});

dem = DEM.Z;
cs  = DEM.cellsize;

% any nans
INAN = isnan(dem);
flaginan = any(INAN(:));

% structuring element
radiuspx = ceil(p.Results.radius/cs);
SE = strel('disk',radiuspx,p.Results.N);


switch p.Results.type
    case 'max'
        % Maximum filter
        if flaginan
            dem(INAN) = -inf;
        end
        H = imdilate(dem,SE);
    case 'min'
        % Minimum filter
        if flaginan
            dem(INAN) = inf;
        end
        H = imerode(dem,SE);
    case 'range'
        if flaginan
            dem(INAN) = -inf;
        end
        H1 = imdilate(dem,SE);
        if flaginan
            dem(INAN) = inf;
        end
        H2 = imerode(dem,SE);
        H  = H1-H2;
    case {'mean','average'}
        if flaginan
            [~,L] = bwdist(~INAN,'e');
            dem = dem(L);
        end            
        H   = fspecial('disk',radiuspx);
        H   = imfilter(dem,H,'symmetric','same','conv');
    case 'median'
        if flaginan
            [~,L] = bwdist(~INAN,'e');
            dem = dem(L);
        end
        H   = getnhood(SE);
        n   = round(sum(H(:))/2);
        H   = ordfilt2(dem,n,H,'symmetric');
        
    case 'prctile'
        if flaginan
            [~,L] = bwdist(~INAN,'e');
            dem = dem(L);
        end
        H   = getnhood(SE);
        n   = round(sum(H(:))*p.Results.prc/100);
        H   = ordfilt2(dem,n,H,'symmetric');
    
    case 'std'
        if flaginan
            [~,L] = bwdist(~INAN,'e');
            dem = dem(L);
        end
        H   = getnhood(SE);
        H   = stdfilt(dem,H);        
end
            

% Handle nans
if flaginan
    H(INAN) = nan;
end

% prepare output
OUT = DEM;
OUT.Z = H;
OUT.name = ['local topography (' p.Results.type ')'];

