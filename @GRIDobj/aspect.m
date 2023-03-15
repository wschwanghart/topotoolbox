function OUT = aspect(DEM,classify)

%ASPECT angle of exposition from a digital elevation model (GRIDobj)
%
% Syntax
%
%     ASP = aspect(DEM)
%     ASP = aspect(DEM,classify)
%
% Description
%
%     aspect returns the slope exposition of each cell in a digital
%     elevation model in degrees. In contrast to the second output of
%     gradient8 which returns the steepest slope direction, aspect 
%     returns the angle as calculated by surfnorm.
%
% Input
%
%     DEM         digital elevation model (class: GRIDobj)
%     classify    false (default) or true. If true, directions are
%                 classified according to the scheme proposed by
%                 Gomez-Plaza et al. (2001)
%
% Output
% 
%     ASP       aspect in degrees (clockwise from top) or classified grid,
%               if classify is set to true
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     ASP = aspect(DEM);
%     imageschs(DEM,ASP)
%
% References
%
%     Gómez-Plaza, A.; Martínez-Mena, M.; Albaladejo, J. & Castillo, V. M.
%     (2001): Factors regulating spatial distribution of soil water content
%     in small semiarid catchments. Journal of Hydrology, 253, 211 - 226.
%
%
% See also: GRIDobj/gradient8, GRIDobj/reclassify
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017

narginchk(1,2)

if nargin == 1
    classify = false;
else
   
end

OUT     = DEM;
OUT.Z   = [];

if classify
    aspedges = (0:45:360)';
    aspclass = [1 3 5 7 8 6 4 2];
end

% Large matrix support. Break calculations in chunks using blockproc
if numel(DEM.Z)>(5001*5001)
    fun    = @(x) aspfun(x);
    blksiz = bestblk(size(DEM.Z),5000);
    OUT.Z  = blockproc(DEM.Z,blksiz,fun,'BorderSize',[1 1]);
else
    OUT.Z = aspfun(DEM.Z);
end

OUT.name = 'aspect';
OUT.zunit = 'degree';




function ASP = aspfun(Z) 
if isstruct(Z)
    Z = Z.data;
end

[Nx,Ny] = surfnorm(Z);
ASP   = cart2pol(Nx,Ny);
ASP   = mod(90+ASP/pi*180,360);
% ASP   = 270 - ASP/pi*180 -90;

if classify
    
    [~,bin] = histc(ASP(:),aspedges);
    bin     = reshape(bin,size(Z));
    ASPc    = zeros(size(Z),'uint8');
    ASPc(bin>0) = aspclass(bin(bin>0));
    ASP     = ASPc;

end
end

end