function DEMr = resample(DEM,target,method,swapzone)

%RESAMPLE change spatial resolution of a GRIDobj
%
% Syntax
%
%     DEMr = resample(DEM,cellsize)
%     DEMr = resample(DEM,GRID)
%     DEMr = resample(...,method)
%     DEMr = resample(DEM,GRID,method,swapzone)
%
% Description
%
%     resample changes the cellsize of a grid. The function uses the MATLAB
%     function imtransform. If an instance of GRIDobj is supplied as
%     second argument, resample interpolates values in DEM to match the
%     spatial reference of GRID.
%
% Input arguments
%
%     DEM       grid object (GRIDobj)
%     cellsize  cellsize of resampled grid
%     GRID      other grid object
%     method    'bicubic', 'bilinear', or 'nearest' 
%     swapzone  true or false. If true and if DEM and GRID have different
%               projected coordinate systems, the function will attempt to
%               reproject and resample the DEM in one step. Note that this
%               requires the mapping toolbox. In case, the DEM is in a
%               geographic coordinate system, please use the function
%               reproject2utm(DEM,GRID).
%
% Output arguments
%
%     DEMr    grid object (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEMr = resample(DEM,100);
%     imagesc(DEMr)
%
%
% See also: griddedInterpolant, imtransform
%        
% Author:  Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 8. August, 2015 

% check input arguments
narginchk(2,4)
validateattributes(target,{'double' 'GRIDobj'},{'scalar'})
if nargin == 2;
    method = 'bilinear';
    swapzone = false;
elseif nargin == 3;
    method = validatestring(method,{'bicubic', 'bilinear', 'nearest' });
    swapzone = false;
else
    method = validatestring(method,{'bicubic', 'bilinear', 'nearest' });
end

% check underlying class
if islogical(DEM.Z)
    method = 'nearest';
end

if swapzone && isa(target,'GRIDobj');
    if ~isequal(DEM.georef.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey,...
            target.georef.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)
        warning('UTM zones differ. Will attempt to match.');
        swapzone = true;
    else
        swapzone = false;
    end
end

% get coordinate vectors
[u,v] = getcoordinates(DEM);

% tform
T = maketform('affine',[1 0 0; 0 1 0; 0 0 1]);

% Fillvalues
if isinteger(DEM.Z)
    fillval = 0;
elseif islogical(DEM.Z)
    fillval = 0;
else
    fillval = nan;
end


if isa(target,'GRIDobj')
    % the target is another GRIDobj
    
    % check whether both grids have the same projection
    if swapzone
        mstructsource = DEM.georef.mstruct;
        mstructtarget = target.georef.mstruct;
        T = maketform('custom', 2, 2, ...
                @FWDTRANS, ...
                @INVTRANS, ...
                []);
    end
     
    
    DEMr    = target;
    [xn,yn] = getcoordinates(DEMr);
    
    DEMr.Z = imtransform(DEM.Z,T,method,...
        'Udata',[u(1) u(end)],'Vdata',[v(1) v(end)],...
        'Xdata',[xn(1) xn(end)],'Ydata',[yn(1) yn(end)],...
        'Size',DEMr.size,...
        'FillValues',fillval);
    DEMr.name = [DEM.name ' (resampled)'];
        
else
    csnew   = target;
    DEMr    = GRIDobj([]);
    [DEMr.Z,xn,yn] = imtransform(DEM.Z,T,method,...
        'Udata',[u(1) u(end)],'Vdata',[v(1) v(end)],...
        'Xdata',[u(1) u(end)],'Ydata',[v(1) v(end)],...
        'XYscale',[csnew csnew],...
        'FillValues',fillval);


    % new referencing matrix
    DEMr.refmat = [0 -csnew;...
                   csnew 0; ...
                   xn(1)-csnew ...
                   yn(1)+csnew];
    % size of the resampled grid           
    DEMr.size    = size(DEMr.Z);
    DEMr.cellsize = csnew;
    DEMr.georef = DEM.georef;
end

DEMr.name    = [DEM.name ' (resampled)'];

if ~isempty(DEMr.georef) && ~isa(target,'GRIDobj');
    DEMr.georef.RefMatrix = DEMr.refmat;
    DEMr.georef.Height = DEMr.size(1);
    DEMr.georef.Width  = DEMr.size(2);
    
    DEMr.georef.SpatialRef = refmatToMapRasterReference(DEMr.refmat, DEMr.size);
    
end


%%

    function x = FWDTRANS(u,~)
        % invtrans first
        [lati,long] = minvtran(mstructsource,u(:,1),u(:,2));
        [x,y] = mfwdtran(mstructtarget,lati,long);
        x = [x y];        
    end

    function u = INVTRANS(x,~)
        [lati,long] = minvtran(mstructtarget,x(:,1),x(:,2));
        [x,y] = mfwdtran(mstructsource,lati,long);
        u = [x y];   
        
    end
end






