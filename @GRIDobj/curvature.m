function C = curvature(DEM,ctype,varargin)

%CURVATURE 8-connected neighborhood curvature of a digital elevation model 
%
% Syntax
%
%     C = curvature(DEM)
%     C = curvature(DEM,type)
%     C = curvature(DEM,type,pn,pv,...)
%
% Description
%     
%     curvature returns the second numerical derivative (curvature) of a
%     digital elevation model. By default, curvature returns the profile
%     curvature (profc). 
%
% Input arguments
%
%     DEM    digital elevation model (GRIDobj)
%     type   'profc' (default) : profile curvature [m^(-1)]
%            'planc' : planform curvature or contour curvature [m^(-1)]
%            'tangc' : tangential curvature [m^(-1)]
%            'meanc' : mean curvature [m^(-1)]
%            'total' : total curvature [m^(-2)]
%
%     Parameter name value/pairs
%     
%     'useblockproc'    true or {false}: use block processing 
%                       (see function blockproc)
%     'useparallel'     true or {false}: use parallel computing toolbox
%     'blocksize'       blocksize for blockproc (default: 5000)
%     'meanfilt'        true or {false}: if true, preprocess DEM with 
%                       [3x3] mean filter. 
%
% Output arguments
%
%     C      curvature (GRIDobj)
%
% Remarks
%     
%     Please note that curvature is not defined for cells with zero 
%     gradient. Here, curvature is set to zero.
%
%     All formulas are according to Schmidt et al. (2003) on page 800.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEM = filter(DEM);
%     C = curvature(DEM,'planc');
%     imageschs(DEM,C,'percentclip',0.1)
%
% Reference
%
%     Schmidt, J., Evans, I.S., Brinkmann, J., 2003. Comparison of
%     polynomial models for land surface curvature calculation.
%     International Journal of Geographical Information Science 17,
%     797-814. doi:10.1080/13658810310001596058
%
% See also: GRIDobj/gradient8
%        
% Author:  Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017


% check input arguments
narginchk(1,inf);
if nargin == 1
    ctype = 'profc';
else
    ctype = validatestring(ctype,{'profc','planc','tangc','meanc','total'});
end

p = inputParser;
p.FunctionName = 'GRIDobj/curvature';
addParamValue(p,'useblockproc',false,@(x) isscalar(x));
addParamValue(p,'blocksize',5000,@(x) isscalar(x));
addParamValue(p,'useparallel',false,@(x) isscalar(x));
addParamValue(p,'meanfilt',false,@(x) isscalar(x));
parse(p,varargin{:});

if p.Results.meanfilt
    DEM = filter(DEM);
end

% create a copy of the DEM instance
C = DEM;
c = class(DEM.Z);
switch c
    case 'double'
        C.Z = double.empty(0,0);
    otherwise
        C.Z = single.empty(0,0);
end

% Large matrix support. Break calculations in chunks using blockproc
% Parallisation for large grids using blockproc does in my experience with
% four cores hardly increase the speed. 
if p.Results.useblockproc
    blksiz = bestblk(size(DEM.Z),p.Results.blocksize); 
    cs  = C.cellsize;
    fun = @(x) curvaturesub(x,cs,ctype); 
    C.Z = blockproc(DEM.Z,blksiz,fun,...
           'BorderSize',[1 1],...
           'Padmethod','symmetric',...
           'UseParallel',p.Results.useparallel);
else
    C.Z = curvaturesub(DEM.Z,C.cellsize,ctype);
end

C.name = ctype;

end
% subfunction

function curv = curvaturesub(dem,cs,ctype)

if isstruct(dem);
    dem = dem.data;
    % DEM has already been padded
    correctedges = false;
    shape = 'same';
else
    correctedges = true;
    shape = 'valid';
end

if correctedges
    dem = padarray(dem,[1 1],'symmetric');
end

% First-order partial derivatives:
% kernel for dz/dx
kernel = [-1 0 1; -1 0 1; -1 0 1]./(6*cs);
fx = conv2(dem,kernel,shape);
% kernel for dz/dy
kernel = [1 1 1; 0 0 0; -1 -1 -1]./(6*cs);
fy = conv2(dem,kernel,shape);

% Second order derivatives according to Evans method (see Olaya 2009)
%
% z1 z2 z3
% z4 z5 z6
% z7 z8 z9

% kernel for d2z/dx2
kernel = [1 -2 1; 1 -2 1; 1 -2 1]./(3*cs.^2);
fxx = conv2(dem,kernel,shape);
% kernel for d2z/dy2
kernel = kernel';
fyy = conv2(dem,kernel,shape);
% kernel for d2z/dxy
kernel = [-1 0 1; 0 0 0; 1 0 -1]./(4*cs.^2);
fxy = conv2(dem,kernel,shape);


%% Other options to calculate Second-order partial derivatives:
% r = gradient(p,cs);
% t = gradient(q',cs)';
% % Second-order mixed partial derivative:
% s = gradient(p',cs)';


switch ctype
    case 'profc'
        curv = - (fx.^2 .* fxx + 2*fx.*fy.*fxy + fy.^2.*fyy)./((fx.^2 + fy.^2).*(1 + fx.^2 + fy.^2).^(3/2));
    case 'tangc'
        curv = - (fy.^2 .* fxx - 2*fx.*fy.*fxy + fx.^2.*fyy)./((fx.^2 + fy.^2).*(1 + fx.^2 + fy.^2).^(1/2));
    case 'planc'
        curv = - (fy.^2 .* fxx - 2*fx.*fy.*fxy + fx.^2.*fyy)./((fx.^2 + fy.^2).^(3/2));
    case 'meanc'
        curv = - ((1+fy.^2).*fxx - 2.*fxy.*fx.*fy + (1+fx.^2).*fyy)./ ...
            (2.* (fx.^2+fy.^2+1).^(3/2));
%         curv = (fx.^2 .* fxx + 2*fx.*fy.*fxy + fy.^2.*fyy)./((fx.^2 + fy.^2).*(1 + fx.^2 + fy.^2)) ...
%             - ((1+fy).^2 .* fxx + 2*fx.*fy.*fxy + (1+fx).^2.*fyy)./(2.*(1 + fx.^2 + fy.^2).^(3/2));
    case 'total'
        curv = fxx.^2 + 2*fxy.^2+fyy.^2;
end

if correctedges
    dem = dem(2:end-1,2:end-1);
end
curv(isinf(curv) | isnan(curv)) = 0;
curv(isnan(dem)) = nan;
curv = reshape(curv,size(dem));

end