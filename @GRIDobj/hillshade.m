function OUT2 = hillshade(DEM,varargin)

%HILLSHADE create hillshading from a digital elevation model (GRIDobj)
%
% Syntax
%    
%     H = hillshade(DEM)
%     H = hillshade(DEM,'pn','pv',...)
%
% Description
%
%     Hillshading is a very powerful tool for relief depiction.
%     hillshade calculates a shaded relief for a digital elevation model 
%     based on the angle between the surface and the incoming light beams.
%     If no output arguments are defined, the hillshade matrix will be
%     plotted with a gray colormap. The hillshading algorithm follows the
%     logarithmic approach to shaded relief representation of Katzil and
%     Doytsher (2003).
%
% Input
%
%     DEM       Digital elevation model (class: GRIDobj)
%
% Parameter name/value pairs
%
%     'azimuth'         azimuth angle, (default=315)
%     'altitude'        altitude angle, (default=60)
%     'exaggerate'      elevation exaggeration (default=1). Increase to
%                       pronounce elevation differences in flat terrain
%     'useblockproc'    true or {false}: use block processing 
%                       (see function blockproc)
%     'useparallel'     true or {false}: use parallel computing toolbox
%     'blocksize'       blocksize for blockproc (default: 5000)
%     'method'          'surfnorm' (default) or 'mdow'
%
%
% Output
%
%     H         shaded relief (ranges between 0 and 1)
%
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     hillshade(DEM)
% 
% References
%
%     Katzil, Y., Doytsher, Y. (2003): A logarithmic and sub-pixel approach
%     to shaded relief representation. Computers & Geosciences, 29,
%     1137-1142.
%
% See also: SURFNORM, IMAGESCHS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. August, 2017



% Parse inputs
p = inputParser;
p.StructExpand  = true;
p.KeepUnmatched = false;
p.FunctionName = 'hillshade'; 
addParamValue(p,'azimuth',315,@(x) isscalar(x) && x>= 0 && x<=360);
addParamValue(p,'altitude',60,@(x) isscalar(x) && x>= 0 && x<=90);
addParamValue(p,'exaggerate',1,@(x) isscalar(x) && x>0);
addParamValue(p,'useparallel',true);
addParamValue(p,'blocksize',2000);
addParamValue(p,'useblockproc',true,@(x) isscalar(x));
addParamValue(p,'method','default');
parse(p,varargin{:});

OUT     = DEM;
OUT.Z   = [];

cs      = DEM.cellsize;
azimuth = p.Results.azimuth;
altitude = p.Results.altitude;
exaggerate = p.Results.exaggerate;
method   = validatestring(p.Results.method,{'default','surfnorm','mdow'});

% Large matrix support. Break calculations in chunks using blockproc
if numel(DEM.Z)>(10001*10001) && p.Results.useblockproc
    blksiz = bestblk(size(DEM.Z),p.Results.blocksize);    
    padval = 'symmetric';
    Z      = DEM.Z;
    % The anonymous function must be defined as a variable: see bug 1157095
    fun   = @(x) hsfun(x,cs,azimuth,altitude,exaggerate,method);
    HS = blockproc(Z,blksiz,fun,...
                'BorderSize',[1 1],...
                'padmethod',padval,...
                'UseParallel',p.Results.useparallel);
    OUT.Z = HS;
else
    OUT.Z = hsfun(DEM.Z,cs,azimuth,altitude,exaggerate,method);
end

OUT.name = 'hillshade';
OUT.zunit = '';

if nargout == 0
    OUT.Z = uint8(OUT.Z*255);
    imagesc(OUT);
    colormap(gray)
else
    OUT2 = OUT;
end

end
%% Subfunction
function H = hsfun(Z,cs,azimuth,altitude,exaggerate,method)

if isstruct(Z)
    Z = Z.data;    
end

switch method
    case {'default','surfnorm'}

        % correct azimuth so that angles go clockwise from top
        azid = azimuth-90;
        
        % use radians
        altsource = altitude/180*pi;
        azisource = azid/180*pi;
        
        % calculate solar vector
        [sx,sy,sz] = sph2cart(azisource,altsource,1);
        
        % calculate surface normals
        [Nx,Ny,Nz] = surfnorm(Z/cs*exaggerate);
        
        % calculate cos(angle)
        % H = [Nx(:) Ny(:) Nz(:)]*[sx;sy;sz];
        % % reshape
        % H = reshape(H,size(Nx));
        
        H = Nx*sx + Ny*sy + Nz*sz;
        
        % % usual GIS approach
        % H = acos(H);
        % % force H to range between 0 and 1
        % H = H-min(H(:));
        % H = H/max(H(:));

    case 'mdow'

        
        % correct azimuth so that angles go clockwise from top
        azid = [360 315 225 270] - 90;
        azisource = azid/180*pi;
        
        altsource = 30/180*pi;
        altsource = repmat(altsource,size(azisource));
        
        % calculate solar vector
        [sx,sy,sz] = sph2cart(azisource,altsource,1);
        
        % calculate surface normals
        [Nx,Ny,Nz] = surfnorm(Z/cs*exaggerate);
        
        
        H = bsxfun(@times,Nx(:),sx) + bsxfun(@times,Ny(:),sy) + bsxfun(@times,Nz(:),sz);
        H = max(H,0);
        H = sum(H,2)./3;
        H = reshape(H,size(Z));
end


end
