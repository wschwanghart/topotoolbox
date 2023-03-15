function varargout = contour(DEM,varargin)

%CONTOUR contour plot of an instance of GRIDobj
%
% Syntax
%
%     contour(DEM,n)
%     contour(DEM,levels)
%     MS = ...
%     [x,y,z] = ...
%
% Description
%
%     Overloaded method for the builtin function contour. If called with
%     output arguments, the function will not plot the contours.
%
% Input arguments
%
%     DEM       grid (class: GRIDobj)
%     n         number of contour levels
%     levels    vector with levels (if only one specific level should be 
%               returned, use [level level]).
%  
% Output arguments
%
%     MS        structure array (can be exported using shapewrite or
%               plotted using mapshow, requires the Mapping Toolbox)
%     [x,y,z]   nan-punctuated vectors for 2d and 3d plotting
%
% Example
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     imageschs(DEM,gradient8(DEM));
%     [x,y] = contour(DEM,10);
%     hold on
%     plot(x,y,'k')
%     hold off
%
%
% See also: GRIDobj, GRIDobj/imageschs, GRIDobj/griddedcontour, contourc
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 6. August, 2013

[Z,x,y] = GRIDobj2mat(DEM);
if nargout == 0   
    contour(x,y,double(Z),varargin{:});
    return
end

c = contourc(double(x),double(y),double(Z),double(varargin{:}));
IXs = 2;
counter = 1;
while IXs(counter) < size(c,2)
    elev(counter)   = c(1,IXs(counter)-1);
    IXe(counter) = IXs(counter) + c(2,IXs(counter)-1) - 1;   
    counter = counter+1;
    IXs(counter) = IXe(counter-1)+2;
    
end
IXs(end) = [];

if isempty(IXs);
    if nargout == 1;
        varargout{1} = struct('Geometry','Line',...
        'X',{},...
        'Y',{},...
        'elev',{});
    else
        varargout{1} = [];
        varargout{2} = [];
        varargout{3} = [];
    end
    return
end
    

nrcontours = numel(IXe);
xy = [];
for r = 1:nrcontours
    
    cc = c(:,IXs(r):IXe(r))';
    xy = [xy; cc; [nan nan]];
end

x = xy(:,1);
y = xy(:,2);

% prepare output
if nargout >= 2 
    
    varargout{1} = x;
    varargout{2} = y;
    
    if nargout == 3
        IXs = IXs-1;
        IXe = IXe-1;
        z = nan(size(x));
        for r = 1:nrcontours
            z(IXs(r):IXe(r)) = elev(r);
        end
        varargout{3} = z;
    end
    
elseif nargout == 1
    
    IXs = num2cell(IXs-1);
    IXe = num2cell(IXe-1);
    
    varargout{1} = struct('Geometry','Line',...
        'X',cellfun(@(s,e) x(s:e),IXs,IXe,'UniformOutput',false),...
        'Y',cellfun(@(s,e) y(s:e),IXs,IXe,'UniformOutput',false),...
        'elev',num2cell(elev));
end
    



