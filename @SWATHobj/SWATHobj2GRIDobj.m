function [varargout] = SWATHobj2GRIDobj(SW,DEM,varargin)
%SWATHOBJ2GRIDOBJ create a GRIDobj with swath-specific information
%
% Syntax
%
%     OUT = SWATHobj2GRIDobj(SW,DEM)
%     OUT = SWATHobj2GRIDobj(SW,DEM,type)
%     [OUT1,OUT2,...] = SWATHobj2GRIDobj(SW,DEM,type1,type2,...)
%
% Description
%
%     SWATHobj2GRIDobj(SW,DEM) creates a GRIDobj the same size as DEM with
%     pixels labeled 1 that are within the swath represented by SWATHobj
%     and all other pixels labeled 0.
%     SWATHobj2GRIDobj(SW,DEM,type) with type being a string that specifies
%     the type of output. type = 'ix' yields the same output as calling
%     the function without any type (see above). type = 'distx' writes to
%     the pixels the along-swath distance. type = 'disty' writes to the
%     pixels the across-swath distance. type = 'z' writes down the z-values
%     of the SWATHobj. If the function is called with several type strings, 
%     several outputs (GRIDobj's) will be created.
%
% Input arguments
%
%     SW     instance of SWATHobj
%     DEM    instance of GRIDobj
%     type   string that specifies the type of output. Valid choices are
%            'distx', 'disty', 'ix', and 'z'. 'distx' and 'disty' give the 
%            distance along and across the swath, respectively. 'ix' 
%            creates a binary image indicating pixels inside the SWATHobj.
%            'z' writes down the z-values of the SWATHobj.
%
% Output arguments
%
%     OUT    instance of GRIDobj
%
% Examples
%
%     Example 1: Label pixels inside SWATHobj
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     x = [3.8272e5, 3.9058e5, 4.0175e5, 4.0680e5];
%     y = [3.7955e6, 3.8017e6, 3.8002e6, 3.7924e6];
%     SW = SWATHobj(DEM,x,y,'smooth',1e3,'width',3e3);
%     IX = SWATHobj2GRIDobj(SW,DEM);
%     figure,imageschs(DEM,IX)
% 
%     Example 2: Label pixels by distance along river 
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     A  = flowacc(FD);
%     S = STREAMobj(FD,A>1000);
%     S = trunk(klargestconncomps(S,1));
%     [x,y] = STREAMobj2XY(S);
%     ix = ~isnan(x);
%     xval = flipud(x(ix));
%     yval = flipud(y(ix));
%     SW = SWATHobj(DEM,xval,yval,'smooth',1e3,'width',3e3);
%     SW = tidy(SW);
%     DX = SWATHobj2GRIDobj(SW,DEM,'distx');
%     figure,imageschs(DEM,DX), hold on
%     plot(SW,'points',false,'labeldist',1e3), hold off
%     colorbar
%
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: May, 2015; updated August, 2017


narginchk(2,5)
if nargin==2 && nargout>1
    error('Too many outputs.')
elseif nargin>2 && length(varargin)~=nargout
    error('Number of inputs mismatches number of outputs.')
end

if nargin == 2
    type{1} = 'ix';
else
    for i = 1 : length(varargin)
        validatestring(varargin{i},{'distx','disty','ix','z'});
    end
    type = varargin;
end


[x,y] = getcoordinates(DEM);
[X,Y] = meshgrid(x,y);
maxx = max(x);
minx = min(x);
maxy = max(y);
miny = min(y);

% Delineate SWATHobj in GRID domain
IM = true(size(SW.Z));
IM(isnan(SW.Z)) = false;
B = bwboundaries(IM,4);
A = false(DEM.size);
for k = 1 : length(B)
    ix_bdy = sub2ind(size(IM),B{k}(:,1),B{k}(:,2));
    x_bdy = SW.X(ix_bdy);
    y_bdy = SW.Y(ix_bdy);
    if length(x_bdy)>2
        [newx,newy] = eqdistline(x_bdy,y_bdy,DEM.cellsize/3);
    else
        newx = x_bdy; newy = y_bdy;
    end
    % take points outside the DEM into account
    newx(newx>maxx) = maxx;
    newx(newx<minx) = minx;
    newy(newy>maxy) = maxy;
    newy(newy<miny) = miny;
    % convert to linear indices
    ix_edge = coord2ind(DEM,newx,newy);
    ix_edge = unique(ix_edge);
    ix_edge = ix_edge(~isnan(ix_edge));
    A(ix_edge) = true;
end
A = imfill(A,4,'holes');

% Swath locations
ix_SW = ~isnan(SW.Z);
x = SW.X(ix_SW);
y = SW.Y(ix_SW);

% Grid locations
ix_grid = find(A==1);
xq = X(ix_grid);
yq = Y(ix_grid);
[xq_edges,yq_edges] = ind2coord(DEM,ix_edge);

varargout = cell(size(type));
for i = 1 : length(type)
    if strcmp(type{i},'ix')
        OUT = DEM;
        OUT.Z(:) = false;
        OUT.Z(ix_grid) = true;
        OUT.Z(ix_edge) = true;
        varargout{i} = OUT;
    else
        switch type{i}
            case 'distx' % longitudinal distance
                D = repmat(SW.distx',length(SW.disty),1);
            case 'disty' % transverse distance
                D = repmat(SW.disty,1,length(SW.distx));
            case 'z' % z-values of the SWATHobj
                D = SW.Z;
        end
        OUT = DEM;
        OUT.Z(:) = nan;
        % interpolate between given points
        v = D(ix_SW);
        F = TriScatteredInterp(x,y,v);
        vq = F(xq,yq);
        OUT.Z(ix_grid) = vq;
        % the edges might be a bit problematic
        F = TriScatteredInterp(x,y,v,'nearest');
        vq = F(xq_edges,yq_edges);
        OUT.Z(ix_edge) = vq;
        varargout{i} = OUT;
    end
end

