function P = polygon2GRIDobj(DEM,MS,field)

%POLYGON2GRIDobj convert polygon to a grid
%
% Syntax
%
%     P = polygon2GRIDobj(DEM,MS)
%     P = polygon2GRIDobj(DEM,MS,field)
%
% Description
%
%     polygon2GRIDobj maps polygons in the mapping structure MS to a 
%     GRIDobj with the same extent and resolution as the GRIDobj DEM.
%
%     Note that polygon2GRIDobj uses polyshape methods that have been
%     released with MATLAB 2017b. The function runs for older versions,
%     too, but may return different results if polygons in MS extend beyond
%     the boundaries of DEM and contain holes.
%     
% Input arguments
%
%     DEM    grid 
%     MS     mapstruct of a polyline as imported by shaperead. Must have
%            the fields X and Y.
%     field  field name of the numeric data to be mapped. If not
%            provided, polygon2GRIDobj will map logical values.
%
% Output arguments
%
%     P      grid (GRIDobj). Grid has the same extent and cellsize
%            as DEM.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','carve');
%     D   = drainagebasins(FD);
%     MS  = GRIDobj2polygon(D);
%     P   = polygon2GRIDobj(D,MS,'ID');
%
%
% Note: This function has not yet been fully tested. Please report bugs.
%
%
% See also: GRIDobj/coord2ind, GRIDobj/sub2coord, GRIDobj/getcoordinates,
%           GRIDobj/createmask, line2GRIDobj, GRIDobj2polygon
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. May, 2018


% 2 or 3 input arguments?
if nargin == 2
    P = GRIDobj(DEM,'logical');
    writelogical = true;
    writeclass = @logical;
else
    % check the class of the field
    cl = class(MS(1).(field));
    P  = GRIDobj(DEM,cl);
    writelogical = false;
    switch cl
        case {'double','single'}
            P = P*nan(1,cl);
        case 'logical'
            writelogical = true;
    end 
    writeclass = @(x) cast(x,cl);
end

% get DEM coordinates
[X,Y] = getcoordinates(DEM);
siz   = DEM.size;

IX_hole = [];

%which matlab version do we have. If MATLAB 9,3
if ~verLessThan('matlab','9.3')
    [xx,yy]  = getoutline(DEM);
    poutline = polyshape(xx,yy);
    poutline = polybuffer(poutline,DEM.cellsize/4);
    warning off
    for r = 1:numel(MS)
        
        psh = polyshape(MS(r).X,MS(r).Y);
        psh = intersect(psh,poutline);
        
        if isempty(psh.Vertices)
            % The feature is completely outside the DEM boundaries
            MS(r).X = [];
            MS(r).Y = [];
            continue
        end
        
        pshhole = holes(psh);
        psh = rmholes(psh);
        % Are the polygons multipart?
        psh = regions(psh);
        
        % loop through regions
        for r2 = 1:numel(psh)
            if r2 == 1
                nrnew = r;
            else
                nrnew = numel(MS);
                nrnew = nrnew + 1;
            end
            xy  = psh(r2).Vertices;
            MS(nrnew).X = xy(:,1);
            MS(nrnew).Y = xy(:,2);
            
            if ~writelogical
                MS(nrnew).(field) = MS(r).(field);
            end
        end
            
        % then loop through holes
        for r2 = 1:numel(pshhole)
            xy  = pshhole(r2).Vertices;
            nrnew = numel(MS);
            nrnew = nrnew + 1;
            MS(nrnew).X = xy(:,1);
            MS(nrnew).Y = xy(:,2);
            IX_hole = [IX_hole nrnew]; %#ok<AGROW>

        end  
    end
    warning on
else
    IX_hole = [];
end

is_hole = false(size(MS));
is_hole(IX_hole) = true;

% loop through features of mapping structure MS
h = waitbar(0);
for r = 1:numel(MS)
    waitbar(r/numel(MS),h,...
        ['Please wait (' num2str(r) '/' num2str(numel(MS)) ')']);
    
    % get coordinates
    x = MS(r).X;
    y = MS(r).Y;
    
    if isempty(x)
        continue
    end
    
    I = isnan(x) | isnan(y);
    x(I) = [];
    y(I) = [];
  
    
    % if the value of that attribute should be written to the grid, this
    % value is extracted here.
    if ~writelogical
        
        if ~is_hole(r)
            val = single([MS(r).(field)]);
        else
            val = single(0);
        end
    else
        if is_hole(r) 
            val = false;
        else
            val = true;
        end
        
    end
    
    % convert coordinates to rows and columns
    [row,col] = coord2pixel(x,y,X,Y);
    % and extract a subset of the image and apply poly2mask to that subset.
    % That is much faster than calling poly2mask for the whole image
    [BW,ext]    = getmask(row,col,siz);
    % Then, write that data back to the main grid P
    FillMat     = P.Z(ext(1):ext(2),ext(3):ext(4)); 
    FillMat(BW) = writeclass(val);
    
    P.Z(ext(1):ext(2),ext(3):ext(4)) = FillMat;

end
close(h)
end

function [BW,ext] = getmask(r,c,siz)
ext = getextent(r,c);
ext = [max(floor(ext(1)),1) min(ceil(ext(2)),siz(1)) ...
       max(floor(ext(3)),1) min(ceil(ext(4)),siz(2))];
sizcrop = [ext(2)-ext(1)+1   ext(4)-ext(3)+1];
rcrop   = r - ext(1) + 1;
rcrop   = max(rcrop,1);
rcrop   = min(rcrop,sizcrop(1));

ccrop   = c - ext(3) + 1;
ccrop   = max(ccrop,1);
ccrop   = min(ccrop,sizcrop(2));

BW = poly2mask(ccrop,rcrop,sizcrop(1),sizcrop(2));
end
        
        
       

function [r,c] = coord2pixel(x,y,X,Y)

% force column vectors
x = x(:);
y = y(:);

dx  = X(2)-X(1);
dy  = Y(2)-Y(1);

c = (x-X(1))./dx + 1;
r = (y-Y(1))./dy + 1;
end


function ext = getextent(x,y)

ext = [min(x) max(x) min(y) max(y)];
end





