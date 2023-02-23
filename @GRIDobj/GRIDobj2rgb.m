function [RGB,x,y] = GRIDobj2rgb(DEM,varargin)

%GRIDobj2RGB Convert GRIDobj to RGB image
%
% Syntax
%
%     [RGB,x,y] = GRIDobj2rgb(DEM)
%     ... = GRIDobj2rgb(DEM,'pn',pv,...)
%
% Description
%
%     GRIDobj2rgb converts a GRIDobj to an RGB image.
%
% Input arguments
%
%     DEM    GRIDobj
%
%     Parameter name/value pairs
%
%     'colormap'     Colormap. Default is parula.
%     'nancolor'     Color to be used for nans
%     'clim'         Color range 
%
% Output arguments
%
%     RGB    n x m x 3 array with uint8 RGB values between 0 and 255 
%     x,y    coordinate vectors
%
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     D = drainagebasins(FD,S);
%     DEM = clip(DEM,D>0);
%     [RGB,x,y] = GRIDobj2rgb(DEM,'colormap',landcolor,...
%                                 'nancolor',[0 0.4470 0.7410]);
%     image(x,y,RGB)
%     axis xy
%
% Example 2: Fuse image with hillshade
%
%     RGB2 = imageschs(DEM,[],'colormap',[1 1 1]);
%     RGB3 = uint8(single(RGB) .* single(RGB2)/255);
%     image(x,y,RGB3)
%
% Example 3: Add cast shadows
%
%     C = castshadow(DEM,'azi',135,'alt',15);
%     RGB4 = GRIDobj2rgb(C,'colormap',ttscm("grayC",255,[0 70]));
%     RGB5 = uint8(single(RGB3) .* (single(RGB4)/255));
%     image(x,y,RGB5)
%
% See also: imageschs, GRIDobj2im
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 21. February, 2023

% Parse input
p = inputParser;
p.FunctionName   = 'GRIDobj2rgb';
addParameter(p,'colormap',parula,@(x) ...
    validateattributes(x,'numeric',{'2d','ncols',3,'>=',0,'<=',1}))
addParameter(p,'nancolor',[1 1 1],@(x) ...
    validateattributes(x,'numeric',{'numel',3}))
addParameter(p,'clim',[],@(x) ...
    validateattributes(x,'numeric',{'increasing','numel',2}))
parse(p,varargin{:})

% Any nans?
I = isnan(DEM);

% Colormap
clr = p.Results.colormap;
Z   = DEM.Z;

if isempty(p.Results.clim)
    
    minval = min(Z(:));
    maxval = max(Z(:));
else
    
    minval = p.Results.clim(1);
    maxval = p.Results.clim(2);
end


valuesatclr = linspace(minval,...
                       maxval,...
                       size(clr,1))';

RGB = interp1(valuesatclr,clr,Z(:),'linear'); 

% map nancolor to nan-pixels
RGB(I.Z(:),:) = repmat(p.Results.nancolor(:)',nnz(I.Z),1);

% apply limits
if ~isempty(p.Results.clim)
    I = DEM(:)<minval;
    RGB(I.Z(:),:) = repmat(clr(1,:),nnz(I.Z),1);
    I = DEM(:)>maxval;
    RGB(I.Z(:),:) = repmat(clr(end,:),nnz(I.Z),1);
end

RGB = reshape(RGB,DEM.size(1),DEM.size(2),3);

RGB = RGB*255;
RGB = uint8(RGB);

if nargout > 1
    [x,y] = getcoordinates(DEM);
end
