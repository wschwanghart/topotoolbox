function varargout = evansslope(DEM,varargin)

%EVANSSLOPE Calculate surface slope using Evans method
%
% Syntax
%
%     G = evansslope(DEM)
%     G = evansslope(DEM,'pn',pv,...)
%     [Gx,Gy] = evansslope(DEM)
%
% Description
%
%     Evans method fits a second-order polynomial to 3x3 subgrids. The
%     parameters of the polynomial are the partial derivatives which are
%     used to calculate the surface slope G = sqrt(Gx^2+Gy^2).
%
%     Evans method approximates the surface by regression surfaces.
%     Gradients are thus less susceptible to noise in the DEM.
%
% Input arguments
%
%     DEM    Digital elevation model (GRIDobj)
%     
%     Parameter name/value pairs
%
%     'padval'    'none' (values along edges will have nans), other options
%                  see padarray. Default is 'replicate'.
%     
% Output arguments
%
%     G      Slope (GRIDobj)
%     Gx,Gy  Partial derivatives
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     G = evansslope(DEM);
%     imageschs(DEM,G,'caxis',[0 1])
%
%
% See also: GRIDobj/gradient8, GRIDobj/arcslope, GRIDobj/curvature
%        
% Author:  Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. January, 2023

% Parse input arguments
p = inputParser;
p.FunctionName = 'GRIDobj/evansslope';
addParameter(p,'padval','replicate')
parse(p,varargin{:})

if ischar(p.Results.padval) || isstring(p.Results.padval)
    padval = validatestring(p.Results.padval,{'none','replicate','symmetric','circular'});
else
    padval = p.Results.padval;
end

switch lower(padval)
    case 'none'
        dem = DEM.Z;
    otherwise   
        dem = padarray(DEM.Z,[1 1],padval);
end

cs    = DEM.cellsize;
shape = 'valid';
kernel = [-1 0 1; -1 0 1; -1 0 1]./(6*cs);
fx = conv2(dem,kernel,shape);
% kernel for dz/dy
kernel = [1 1 1; 0 0 0; -1 -1 -1]./(6*cs);
fy = conv2(dem,kernel,shape);

if nargout == 1
    
    switch lower(p.Results.padval)
        case 'none'
            varargout{1} = GRIDobj(DEM)*nan;
            varargout{1}.Z(2:end-1,2:end-1) = sqrt(fx.^2 + fy.^2);
        otherwise
            varargout{1} = GRIDobj(DEM);
            varargout{1}.Z = sqrt(fx.^2 + fy.^2);
    end
    varargout{1}.name = 'Gradient (Evans)';
else
    switch lower(p.Results.padval)
        case 'none'
            varargout{1} = GRIDobj(DEM)*nan;
            varargout{1}.Z(2:end-1,2:end-1) = fx;
            varargout{2} = GRIDobj(DEM)*nan;
            varargout{2}.Z(2:end-1,2:end-1) = fy;
        otherwise
            varargout{1} = DEM;
            varargout{1}.Z = fx;
            varargout{2} = DEM;
            varargout{2}.Z = fy;
    end
    varargout{1}.name = 'fx (Evans)';
    varargout{2}.name = 'fy (Evans)';

end


