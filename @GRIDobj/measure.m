function measure(DEM,varargin)

%MEASURE take interactive measurements along a polyline
%
% Syntax
%
%     measure(DEM)
%     measure(DEM,pn,pv,...)
%
% Description
%
%     This simple gui lets you measure distances and slopes along a
%     polyline. 
%
% Input arguments
%
%     DEM     digital elevation model (class: GRIDobj)
%
%     parameter name, parameter value
%   
%     'xyprecision'   displayed precision of x and y coordinates 
%                     (C-style format string, default '%6.0f')
%     'slopeunit'     unit of slope ({'degree'},'tan','percent') 
%     'showhelp'      show help window ({true} or false)
%     'position'      nx2 matrix with x and y coordinates
%     'colormap'      choose colormap. Default is parula
%     'water'         by default empty. Provide a GRIDobj with water depth
%                     to show water surface and plot and save profiles of
%                     topography and water surface.
% 
% 
% See also: IMDISTLINE, IMROI, IMPOLY, GRIDobj/DEMPROFILE
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 12. February, 2022

% 24. November, 2014: added structure field nodes if profile is saved to
% workspace. This allows to continue working with the measure tool with an
% existing profile by calling 
% measure(DEM,'position',p.nodes)

% 12. February, 2022: add water option to display water surfaces in
% profiles. Increased resolution of profiles


% get DEM coordinates
[X,Y] = getcoordinates(DEM);
X   = X(:);
Y   = Y(:);
dx  = X(2)-X(1);
dy  = Y(2)-Y(1);

xmax  = max(X);
ymax  = max(Y);
xmin  = min(X);
ymin  = min(Y);

% default position
defpos = [xmin + (xmax-xmin)*[1/3; 2/3] ymin + (ymax-ymin)*[2/3; 1/3]];

% choose colormap
if license('test', 'map_toolbox')
    cmap = demcmap([max(DEM) min(DEM)]-1,255);
else
    cmap = parula(255);
end

p = inputParser;

p.FunctionName = 'STREAMobj';
addParamValue(p,'xyprecision','%6.0f',@(x) ischar(x));
addParamValue(p,'slopeunit','degree',@(x) ischar(validatestring(x,{'tan','degree','percent'})));
addParamValue(p,'reset',false);
addParamValue(p,'colormap',cmap);
addParamValue(p,'showhelp',true,@(x) isscalar(x));
addParamValue(p,'position',defpos,@(x) numel(x) >= 4 && size(x,2) == 2 && ...
    (max(x(:,1)) <= xmax && min(x(:,1)) >= xmin && ...
     max(x(:,2)) <= ymax && min(x(:,2)) >= ymin));
addParamValue(p,'water',[]);

parse(p,varargin{:});

if ~p.Results.reset
    if isempty(p.Results.water)
        imageschs(DEM,DEM,'colormap',p.Results.colormap)
    else
        imageschs(DEM,p.Results.water,'colormap',flowcolor)
    end
    zoom reset
end
    
hax = gca;
hf  = gcf;

title('')
fcn = makeConstrainToRectFcn('impoly',get(gca,'XLim'),get(gca,'YLim'));
h = impoly(hax,p.Results.position,'closed',false,'PositionConstraintFcn',fcn);
htext = [];
htextedge = [];
totdist = [];
elev = [];
segdist = [];
switch p.Results.slopeunit
    case 'tan'
        slopefun = @(x) x;
        slopesign = '';
    case 'degree'
        slopefun = @(x) atand(x);
        slopesign = '';
    case 'percent'
        slopefun = @(x) 100*x;
        slopesign = '';
end


textformat = {'HorizontalAlignment','left',...
              'VerticalAlignment','bottom',...
              'BackgroundColor',[0 0 0],...
              'FontSize',8,...
              'Color',[1 1 1]};

textformatedge = {'HorizontalAlignment','center',...
              'VerticalAlignment','middle',...
              'BackgroundColor',[1 1 1],...
              'FontSize',8,...
              'Color',[0 0 0]};   
          
% nr of nodes
showpointtext = false;
nrnodes       = size(getPosition(h),1);          
showedgetext  = nrnodes < 5;

id = addNewPositionCallback(h,@getinfo);
iptcallback = get(hf,'WindowKeyPressFcn');
set(hf,'WindowKeyPressFcn', @toggleview)
% iptcallbackr = get(hf,'WindowKeyReleaseFcn');
% set(hf,'WindowKeyReleaseFcn', @releasefun);

pos = getPosition(h);
getinfo(pos);

if p.Results.showhelp
    showhelp;
end


    function pos = getinfo(pos)
        % GETINFO callback when polyline is modified
        
        nrpos   = size(pos,1);
        
        if nrpos == 1;
            delete(h);
            if ~isempty(htext)
                delete(htext)
            end
            if ~isempty(htextedge)
                delete(htextedge)
            end
            title('measurement tool aborted')
            set(hf,'WindowKeyPressFcn',{})
            return
        end
        
        segdist = sqrt(sum(diff(pos).^2,2));
        totdist = sum(segdist);
        ix      = coord2indsub(pos);
        elev    = DEM.Z(ix);
        slop    = slopefun(-diff(elev)./segdist);
        
        title(['distance: ' num2str(totdist)]);
        
        if showpointtext
            axlim = axis;
            xoffset = axlim(2)-axlim(1);
            yoffset = axlim(4)-axlim(3);
            
            textpos = [pos(:,1)+(xoffset*0.01),pos(:,2)+(yoffset*0.01)];
            
            for r = 1:size(pos,1)
                
                textstr = ['x: ' num2str(pos(r,1),p.Results.xyprecision) ...
                    '\newline' 'y: ' num2str(pos(r,2),p.Results.xyprecision) ...
                    '\newline' 'z: ' num2str(elev(r))];
                try
                    set(htext(r),'Position',textpos(r,:),'String',textstr);
                catch
                    htext(r) = text(textpos(r,1),textpos(r,2),textstr,textformat{:});
                end
            end
            
            if nrpos<numel(htext)
                delete(htext(nrpos+1:end));
                htext(nrpos+1:end) = [];
            end
        else
            delete(htext);
            htext = [];
        end
        
        if showedgetext
            
            for r = 1:size(pos,1)-1
                
                textstredge = ['d: ' num2str(segdist(r),p.Results.xyprecision) ...
                    '\newline' 'slope: ' num2str(slop(r)) slopesign];
                textposedge = (pos(r,:)+pos(r+1,:))/2;
                try
                    set(htextedge(r),'Position',textposedge,'String',textstredge);
                catch
                    htextedge(r) = text(textposedge(1),textposedge(2),textstredge,textformatedge{:});
                end
            end
            if (nrpos-1)<numel(htextedge)
                delete(htextedge(nrpos:end));
                htextedge(nrpos:end) = [];
            end
        else
            delete(htextedge)
            htextedge = [];
        end
        
        
    end


    function ix = coord2indsub(pos)
        % COORD2INDSUB  modified version of GRIDobj/coord2ind
        IX1 = (pos(:,1)-X(1))./dx + 1;
        IX2 = (pos(:,2)-Y(1))./dy + 1;
        
        IX1 = max(min(IX1,DEM.size(2)),1);
        IX2 = max(min(IX2,DEM.size(1)),1);
        
        IX1 = round(IX1);
        IX2 = round(IX2);
        
        
        ix  = sub2ind(DEM.size,IX2,IX1);
    end

    function toggleview(varargin)
        % TOGGLEVIEW callback for key press
        iptcallback(varargin{:});
        switch varargin{2}.Key
            case 't'
                showpointtext = ~showpointtext;
                getinfo(getPosition(h));
            case 'd'
                showedgetext = ~showedgetext;
                getinfo(getPosition(h));
            case 'i'
                zoom(1.5);
            case 'o'
                zoom(2/3);
            case 'r'
                setPosition(h,p.Results.position);
            case 'l'
                zoomtoroi  
                
            case 'p'
                n = ceil((totdist/DEM.cellsize)*5);
                posn = getPosition(h);
                [dn,z] = demprofile(DEM,n,posn(:,1),posn(:,2));
                if ~isempty(p.Results.water)
                    [~,zw] = demprofile(p.Results.water,n,posn(:,1),posn(:,2));
                end
                posfig = get(hf,'OuterPosition');
                
                hfprofile = figure('OuterPosition',posfig.*[1 1 .5 .5]);
                ax = gca;
                plot(ax,dn,z,'k');
                hold on
                if ~isempty(p.Results.water)
                    plot(ax,dn,z+zw,'b');
                end
                plot(ax,cumsum([0;segdist(:)]),elev,'sr');
                hold off
                xlabel('distance [m]');
                ylabel('elevation [m]');
                
            case 's'
                % save profile structure array to workspace
                n = ceil((totdist/DEM.cellsize)*5);
                posn = getPosition(h);
                [S.distance,S.z,S.x,S.y] = demprofile(DEM,n,posn(:,1),posn(:,2));
                if ~isempty(p.Results.water)
                    [~,S.zw] = demprofile(DEM+p.Results.water,n,posn(:,1),posn(:,2));
                end
                S.nodes = posn;
                
                

                answer = inputdlg('Enter variable name:', 'Export profile to workspace', 1, {'profile_struct'});
                try
                if isvarname(answer{1})
                    assignin('base',answer{1},S);
                
                end
                catch
                end

                
                
            case 'h'
                showhelp

        end
    end

%     function releasefun(varargin)
        % TOGGLEVIEW callback for key press
%         iptcallbackr(varargin{:});
%         switch varargin{2}.Key
%             case 'p'
%                 pan off
                
%         end
%     end

    function showhelp
        str = ['Access tools with these keys \n' ...
               'a: add new vertex (move cursor over edge) \n' ...
               't: toggle location information \n' ...
               'd: toggle distance information \n' ...
               'o: zoom out \n' ...
               'i: zoom in \n' ...
               'l: zoom to polyline \n' ...
               'p: plot profile \n' ...
               's: export profile to workspace \n' ...
               'r: reset \n' ...
               'h: show this help'];
         helpdlg(sprintf(str),'measure tool help')
    end
           
    function zoomtoroi
        % zoom to polyline
        posn  = getPosition(h);
        axlim = axis(hax);
        % get x/y axis ratio
        ratio  = (axlim(2)-axlim(1))/(axlim(4)-axlim(3));
        maxpos = max(posn);
        minpos = min(posn);
        ranpos = maxpos-minpos;
        ratiop = ranpos(1)/ranpos(2);
        
        if any(ranpos*1.2 > [xmax-xmin ymax-ymin])
            axis([xmin xmax ymin ymax])
            return
        end
        
        if ratiop >= ratio
            xlims = [minpos(1) maxpos(1)] + [-.1 .1].*ranpos(1); 
            ywidth = (xlims(2)-xlims(1))/ratio;
            ylims = (maxpos(2)+minpos(2))/2 + [-ywidth/2 +ywidth/2];
            if ylims(2) > ymax;
                ylims(2) = ymax;
                ylims(1) = ylims(2)-ywidth;
            end
            if ylims(1) < ymin;
                ylims(1) = ymin;
                ylims(2) = ylims(1)+ywidth;
            end
            
            newaxlim = [xlims ylims];
        else
            ylims = [minpos(2) maxpos(2)] + [-.1 .1].*ranpos(2); 
            xwidth = (ylims(2)-ylims(1))*ratio;
            xlims = (maxpos(1)+minpos(1))/2 + [-xwidth/2 +xwidth/2]; 
            if xlims(2) > xmax;
                xlims(2) = xmax;
                xlims(1) = xlims(2)-xwidth;
            end
            if xlims(1) < xmin;
                xlims(1) = xmin;
                xlims(2) = xlims(1)+xwidth;
            end
            newaxlim = [xlims ylims];
        end
        axis(hax,newaxlim)

    end
end