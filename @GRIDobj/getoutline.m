function varargout = getoutline(DEM,nnan,wm)

%GETOUTLINE get or plot extent of GRIDobj
%
% Syntax
%
%     getoutline(DEM)
%     MS = getoutline(DEM)
%     [x,y] = getoutline(DEM)
%     ... = getoutline(DEM,removenans)
%     getoutline(DEM,removenans,wm)
%
% Description
%
%     getoutline plots the extent of a GRIDobj or returns the coordinate
%     vectors that generate the plot. By default, getoutline returns the
%     coordinate vectors of the DEM edges. By setting removenans = true,
%     you can get the outline around the valid (non-nan) data in the DEM.
%
% Input arguments
%
%     DEM         GRIDobj
%     removenans  set to true, if the outline should be around the valid
%                 data in the DEM.
%     wm          plot in webmap (true or {false})
%
% Output arguments
%
%     MS     mapping structure that can be displayed using mapshow or
%            exported with shapewrite.
%     x,y    coordinate vectors that can be used to plot the extent
%            rectangle (plot(x,y))
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     getoutline(DEM,true,true)
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. October, 2018

nargoutchk(0,2);

if nargin >= 2
    if nnan
        I = isnan(DEM);
        if ~any(I)
            nnan = false;
        end
    end
else
    nnan = false;
end

if ~nnan
    [x,y] = getcoordinates(DEM);
    maxx = max(x);
    minx = min(x);
    
    maxy = max(y);
    miny = min(y);
    
    x = [minx minx maxx maxx minx];
    y = [miny maxy maxy miny miny];
else
    I = ~I;
    B = bwboundaries(I.Z);
    
    x = [];
    y = [];
    for k = 1:length(B)
        boundary = B{k};
        [xb,yb]  = sub2coord(I,boundary(:,1),boundary(:,2));
        x = [x;xb;nan];
        y = [y;yb;nan];
    end
end

% simplify lines
xy = dpsimplify([x(:) y(:)],100*eps);
x = xy(:,1);
y = xy(:,2);


if nargout == 0
    % No output. Plot outline
    if nargin <= 2
        wm = false;
    end
    
    if ~wm
        plot(x,y,'LineWidth',3,'Color',[.8 .8 .8]);
        hh = ishold;
        hold on
        plot(x,y,'k--','LineWidth',1);
        if ~hh
            hold off;
        end
    else
        [lat,lon] = minvtran(DEM.georef.mstruct,x,y);
        wmline(lat,lon)
    end
        
elseif nargout == 1
    % One output argument, create mapping structure
    % split at nans
    if ~isnan(x(end))
        x = [x;nan];
        y = [y;nan];
    end
            
    I = isnan(x);
    nrlines = nnz(I);
    ix = find(I)';
    ixs = [1; ix(1:end-1)'+1];
    ixe = ix-1;
    for r = 1:nrlines
        
        MS(r).Geometry = 'Polygon';
        MS(r).X = x(ixs(r):ixe(r));
        MS(r).Y = y(ixs(r):ixe(r));
        MS(r).ID = r;
        MS(r).name = DEM.name;
    end
    varargout{1} = MS;
elseif nargout == 2
    % Two outputs, return x and y coordinates
    
    varargout{1} = x;
    varargout{2} = y;
end

end
    
       