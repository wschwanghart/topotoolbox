function varargout = plot(D,varargin)
%PLOT    plot the divide network 
%
% Syntax
%
%     plot(D)
%     plot(D,'pn','pv',...)
%     h = plot(...)
%
% Description
%
%     Plot the drainage divide network in map view.
%     
% Input arguments
%
%     D   divides (class DIVIDEobj)
%
%     Parameter name/value pairs
%
%     'style'       'plain' (default if D has no order), 'order' (default 
%                   if D has been ordered) or 'animated'
%     'minwidth'    minimum line width if style is 'order' or
%                   'animated' {1}
%     'maxwidth'    maximum line width if style is 'order' or
%                   'animated' {5}
%     'endpoints'   toggle for plotting the divide network endpoints 
%                   {true} if style is plain, {false} is style is 'order'
%     'junctions'   toggle for plotting the divide network endpoints 
%                   {true} if style is plain, {false} is style is 'order'
%     'pause'       Pause time in seconds if style is 'animated' {0.1}
%     'divlimit'    divide attribute ({'distance','order','none'}) used for
%                   thresholding with values in 'limits'
%     'limit'       min and max values for thresholding the divide network
%                   {1e3 inf}
%
%     Additional parameter name/value pairs for the linear "plot" function 
%     of MATLAB can be provided, although some may not yield the expected
%     outcome if they are in conflict with the above parameters.
%
% Output arguments
%
%     h  graphics handles
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     ST = STREAMobj(FD,flowacc(FD)>1000);
% 
%     figure
%     subplot(3,2,1)
%     imageschs(DEM)
%     title('DEM')
%     colorbar off
%     set(gca,'XTick',[],'YTick',[])
% 
%     subplot(3,2,2)
%     plot(ST,'color','b')
%     axis image
%     title('Drainage network')
%     set(gca,'XTick',[],'YTick',[])
% 
%     subplot(3,2,3)
%     D = DIVIDEobj(FD,ST,'network',false);
%     plot(D)
%     axis image
%     box on
%     title('Drainage divides')
%     set(gca,'XTick',[],'YTick',[])
% 
%     subplot(3,2,4)
%     D = divnet(D,FD);
%     plot(D)
%     axis image
%     box on
%     title('Drainage divide network (DDN)')
%     set(gca,'XTick',[],'YTick',[])
% 
%     subplot(3,2,5)
%     D = sort(D);
%     D = divorder(D);
%     plot(D)
%     axis image
%     box on
%     title('DDN, sorted')
%     set(gca,'XTick',[],'YTick',[])
% 
%     subplot(3,2,6)
%     D = divdist(D);
%     plot(D)
%     axis image
%     box on
%     title('DDN, divide distance > 1000 m')
%     set(gca,'XTick',[],'YTick',[])
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: December 2018

%% Parse inputs
p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'style',[],@(x) ischar(validatestring(x,{'plain','order','animated'})))
addParameter(p,'minwidth',1);
addParameter(p,'maxwidth',5);
addParameter(p,'endpoints',[]);
addParameter(p,'junctions',[]);
addParameter(p,'pause',0.1);
addParameter(p,'divlimit','distance',@(x) ischar(validatestring(x,{'distance','order','none'})));
addParameter(p,'limit',[1e3 inf],@(x) isnumeric(x));

% Parse
parse(p,varargin{:});

% default values
if isempty(p.Results.style) 
    if not(isempty(D.order))
        style = 'order';
    else
        style = 'plain';
    end
else
    style = validatestring(p.Results.style,{'plain','order','distance','animated'});
end

if isempty(p.Results.endpoints) 
    if strcmp(style,'plain')
        endpoints = true;
    else
        endpoints = false;
    end
else
    endpoints = p.Results.endpoints;
end

if isempty(p.Results.junctions) 
    if strcmp(style,'plain')
        junctions = true;
    else
        junctions = false;
    end
else
    junctions = p.Results.junctions;
end

if strcmp(p.Results.divlimit,'none') || (strcmp(p.Results.divlimit,...
        'distance') && isempty(D.distance)) || (strcmp(p.Results.divlimit,...
        'order') && isempty(D.order)) || strcmp(style,'plain')
    ix = true(size(D.IX));    
else
    ix = (D.(p.Results.divlimit)>=min(p.Results.limit) & ...
        D.(p.Results.divlimit)<=max(p.Results.limit)) | isnan(D.IX);
end

if ismember(style,{'order','animated'}) 
    siz = [p.Results.minwidth p.Results.maxwidth];
    do   = D.order(ix);
    uqdo = unique(do(not(isnan(do))));
    maxdo = max(uqdo);
    mindo = min(uqdo);
    a = (siz(2)-siz(1))/(maxdo-mindo);
    b = siz(1)-a*mindo;
end

[x,y] = ind2coord(D,D.IX(ix));

% Create cellarray from p.Unmatched
pn = fieldnames(p.Unmatched);
pv = struct2cell(p.Unmatched);

pnpv = [pn pv];
pnpv = pnpv';
pnpv = pnpv(:)';


axishold = ishold(gca);

switch style
    
    case 'plain'         
        h(1) = plot(x,y,pnpv{:});
        legend(h,'divide');
        
        if junctions && not(isempty(D.jct))
            hold on
            [xjct,yjct] = ind2coord(D,D.jct);
            h(2) = scatter(xjct,yjct,20,'go','filled','DisplayName','junction');
        end
        if endpoints && not(isempty(D.ep))
            hold on
            [xep,yep] = ind2coord(D,D.ep);
            h(3) = scatter(xep,yep,20,'ro','filled','DisplayName','endpoint');
        end
        
        if ~axishold
            hold(gca,'off');
        end
        
    case 'order'
        if isempty(do)
            error('D has no order. Use D = divorder(D) to add orders');
        end
        
        hold on
        for i = length(uqdo):-1:1
            ix = do==uqdo(i);
            h(i) = plot(x(ix),y(ix),pnpv{:},'linewidth',a.*uqdo(i)+b);
        end
        
    case 'animated'
        if isempty(do)
            error('D has no order. Use D = divorder(D) to add orders');
        end
        h = plot(x,y,'color',[.9 .9 .9]);
        hold on
        axis equal
        pause(p.Results.pause);
        for i = length(uqdo):-1:1
            ix = do==uqdo(i);
            h = [h plot(x(ix),y(ix),pnpv{:},'linewidth',a.*uqdo(i)+b)];
            title(['order > ',num2str(i)])
            pause(p.Results.pause);
        end
        hold off
        
end


if nargout==1
    varargout{1} = h;
end

end
