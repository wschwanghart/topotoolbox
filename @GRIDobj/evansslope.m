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
%     'modified'  {false} or true. If true, the surface is weakly smoothed
%                 before gradients are calculated (see Shary et al., 2002)
%     'removenans' [true} or false. Inter- and extrapolate to remove nans 
%                 before calculation. 
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
% Reference: Shary, P. A., Sharaya, L. S., and Mitusov, A. V.: Fundamental 
% quantitative methods of land surface analysis, Geoderma, 107, 1–32, 
% https://doi.org/10.1016/S0016-7061(01)00136-7, 2002.
%
% See also: GRIDobj/gradient8, GRIDobj/arcslope, GRIDobj/curvature
%        
% Author:  Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. January, 2023

% Parse input arguments
p = inputParser;
p.FunctionName = 'GRIDobj/evansslope';
addParameter(p,'padval','replicate')
addParameter(p,'modified',false)
addParameter(p,'removenans',true)
parse(p,varargin{:})

if p.Results.removenans
    I = isnan(DEM);

    if any(I)
        % inpaint
        DEM = inpaintnans(DEM);
        if any(isnan(DEM))
            % extrapolate using nearest neighbor
            [~,L]   = bwdist(~isnan(DEM.Z)); 
            DEM.Z   = DEM.Z(L);
        end
    end
end
      

if p.Results.modified
    kernel = [0 1 0; 1 41 1; 0 1 0]/45;
    DEM.Z = conv2(padarray(DEM.Z,[1 1],'replicate'),kernel,'valid');
end

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
    if exist('I','var')
        varargout{1} = clip(varargout{1},~I);
    end
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
    if exist('I','var')
        I = ~I;
        varargout{1} = clip(varargout{1},I);
        varargout{2} = clip(varargout{2},I);
    end

end


