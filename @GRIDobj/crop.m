function [DEMc,MASK] = crop(DEM,varargin)

%CROP crop an instance of GRIDobj with axis-aligned minimum bounding box
%
% Syntax
%
%     DEMc = crop(DEM)
%     DEMc = crop(DEM,I)
%     DEMc = crop(DEM,I,fillval)
%     DEMc = crop(DEM,ix)
%     DEMc = crop(DEM,x,y)
%     DEMc = crop(DEM,'interactive');
%
% Description 
%
%     crop removes outer parts of a grid based on missing values (NaNs), 
%     a mask, linear indices or coordinate pairs such that the axis-aligned 
%     minimum bounding rectangle remains.
%     
% Input arguments
%
%     DEM      instance of GRIDobj
%     I        mask (logical) as instance of GRIDobj
%     fillval  scalar; values in the output grid where I.Z==0 obtain this 
%              value (crop and clip by setting fillval to nan) 
%     ix       linear index into the DEM
%     x,y      coordinate vectors
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     MASK = createmask(DEM);
%     DEMc = crop(DEM,MASK,NaN);
%     imagesc(DEMc)
%
%
% See also: IND2SUB, GRIDobj/pad, GRIDobj/resample
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 7. March, 2018


narginchk(1,3);

if nargin == 1
    MASK = isnan(DEM.Z);
    
    if ~any(MASK(:))
        DEMc = DEM;
        if nargout == 2
            MASK = DEM;
            MASK.Z = true(DEM.size);
            MASK.zunit = '';
            MASK.name = 'mask for cropping';
        end
        
        return
    end
    MASK = ~MASK;
    % MASK = bwperim(MASK);
    IX  = find(MASK); 
     
elseif nargin >= 2
    if isa(varargin{1},'GRIDobj') || isa(varargin{1},'logical')
        % GRIDobj
        validatealignment(DEM,varargin{1});
        if isa(varargin{1},'GRIDobj')
            MASK = logical(varargin{1}.Z);
        else
            MASK = varargin{1};
        end
        
        if nargin == 3
            validateattributes(varargin{2},{'numeric'},{'scalar'},'crop','fillval',3);
            if isnan(varargin{2}) && ~(isa(DEM.Z,'single') || isa(DEM.Z,'double')) 
                DEM.Z(~MASK) = 0;
                warning('TopoToolbox:GRIDobj',...
                    ['fillval set to zero since data type of the input grid\newline' ...
                     'does not support nans'])
            else
                DEM.Z(~MASK) = varargin{2};
            end
        end
        % MASK = bwperim(MASK);
        IX  = find(MASK);
    elseif nargin == 2
        if isnumeric(varargin{1})
            % indices are supplied
            IX  = varargin{1};
            if numel(IX) < 2
                error('TopoToolbox:GRIDobj',...
                    'At least two indices are required to crop the grid.')
            end
            
            if any(IX<1) || any(IX>prod(DEM.size))
                error('TopoToolbox:GRIDobj',...
                    ['Index must range between 1 and ' num2str(prod(DEM.size)) '.'])
            end
            
        else
            % interactive part
            try
            imagesc(DEM);
            [xx,yy] = getcoordinates(DEM);
            minx = min(xx) - 0.1*DEM.cellsize;
            maxx = max(xx) + 0.1*DEM.cellsize;
            miny = min(yy) - 0.1*DEM.cellsize;
            maxy = max(yy) + 0.1*DEM.cellsize;
            
            c = uicontrol('Style','Text','Units','normalized','position',[0 0 1 0.05],...
                          'String','Draw rectangle and double click when finished.');
            
            fcn = makeConstrainToRectFcn('imrect',[minx maxx],[miny maxy]);
            h = imrect(gca,'PositionConstraintFcn',fcn);
            
            
            addNewPositionCallback(h,@(pos) set(c,'String',...
                ['LX:' num2str(round(pos(1)),'%d') ', LY:' num2str(round(pos(2)),'%d') ...
                 ', UX:' num2str(round(pos(1)+pos(3)),'%d') ', UY:' num2str(round(pos(2)+pos(4)),'%d')]));
            [~] = wait(h);
            MASK  = createMask(h);
            delete(h)
            close
            catch ME
                delete(h);
                error('TopoToolbox:crop','output variable undefined')
                return
            end
            
            % MASK  = bwperim(MASK);
            IX    = find(MASK);
        end
            
    else

        [X,Y] = refmat2XY(DEM.refmat,DEM.size);
        x = varargin{1};
        y = varargin{2};
        
        x = max(min(X),x);
        x = min(max(X),x);
        
        y = max(min(Y),y);
        y = min(max(Y),y);
        
        IX = coord2ind(X,Y,x,y);
    end
    
end


if nargout == 2
    MASK = DEM;
    MASK.Z = false(DEM.size);
    MASK.Z(IX) = true;
    MASK.name = 'mask for cropping';
    MASK.zunit = '';
end


% nr of dimensions
siz  = DEM.size;
k    = [1 cumprod(siz(1:end-1))];

% preallocate subsref structure
S    = substruct('()',cell(1,2));

% subset size
sizout = zeros(1,2);

% loop through dimensions (see ind2sub) 
% and get subscripts of minimum bounding rectangle/box/...
for r = 2:-1:1  
    IX2       = rem(IX-1,k(r))+1;         
    subdim    = (IX-IX2)/k(r)+1; 
    S.subs{r} = min(subdim):max(subdim);
    sizout(r) = numel(S.subs{r});
    IX        = IX2;   
end

DEMc      = DEM;
DEMc.Z    = reshape(subsref(DEM.Z,S),sizout);
DEMc.size = sizout;
[x,y]     = sub2coord(DEM,S.subs{1}(1),S.subs{2}(1));
DEMc.refmat(3,:) = [x y] - DEMc.cellsize*[1 -1];

DEMc.name   = [DEM.name ' (cropped)'];
DEMc.zunit  = DEM.zunit;
DEMc.xyunit = DEM.xyunit;

if ~isempty(DEM.georef)
    % Copy all referencing information
    DEMc.georef = DEM.georef;
    DEMc.georef.SpatialRef = refmatToMapRasterReference(DEMc.refmat, DEMc.size);
    
end
    
    





