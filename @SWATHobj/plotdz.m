function h = plotdz(SW,varargin)
%PLOTDZ creates distance-elevation plot of SWATHobj
%
% Syntax
%
%     plotdz(SW)
%     plotdz(SW,'pn','pv',...)
%     h = ...;
%
% Description
%
%     PLOTDZ creates a profile-view plot of a SWATHobj, showing the
%     statictics (min, max, mean +/- standard dev.) of the z-values as a
%     function of distance along the profile of the SWATHobj.
%
% Input arguments
%
%     SW     instance of SWATHobj
%
%     Parameter name/value pairs
%
%     'left'   {true}, false
%     determines if the left half (as seen from the direction of the 
%     central line) of the SWATHobj is included in the statistics
%
%     'right'   {true}, false
%     determines if the right half (as seen from the direction of the 
%     central line) of the SWATHobj is included in the statistics
%
%     'distadjust'   {0}, scalar
%     allows shifting the x-axis by a scalar value. This is useful when
%     alligning a SWATHobj that was obtained along, e.g., a reach of
%     a drainage network, with the STREAMobj that it was derived from
%
%     'plotminmax'    {true}, false
%     determines whether min/max ranges are plotted
%  
%     'boundedline'   true, {false}
%     both standard deviation from the mean and min/max are plotted using 
%     patch.
%
%     if boundedline is true, then following parameters apply
%     'facecolor'          color (three element vector) for the face color
%                          of the std range. Default is 
%                          [0.5059 0.8471 0.8157]
%     'minmaxfacecolor'    color of the min/max range. Default is [.8 .8 .8]
%     'edgecolor           edge color for the std range {'none'}
%     'minmaxedgecolor'    edge color for the min/max range {'none'}
%     'facealpha'          alpha of the std range {1}
%     'minmaxfacealpha'    alpha of the min/max range {1}
%
% Output arguments (optional)
%
%     h    handle object of the plotted features. h contains several
%          handles to the different parts of the plot
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     SW = SWATHobj(DEM,'dx',200,'dy',200);
%     G = gradient8(DEM,'degree');
%     SWG = mapswath(SW,G);
%     figure, plotdz(SWG)
%     title('Slope angle along swath')
%
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
%         Wolfgang Schwanghart 
% Date: 11. December, 2022

% default patch color
clr = [0.5059 0.8471 0.8157];

% Parse inputs
p = inputParser;
p.FunctionName = 'plotdz';
addRequired(p,'SW',@(x) isa(x,'SWATHobj'));
addParameter(p,'left',true,@(x) islogical(x))
addParameter(p,'right',true,@(x) islogical(x))
addParameter(p,'distadjust',0,@(x) isnumeric(x))
addParameter(p,'boundedline',false,@(x) islogical(x))
addParameter(p,'facecolor',clr)
addParameter(p,'facealpha',1)
addParameter(p,'minmaxfacecolor',[.8 .8 .8])
addParameter(p,'edgecolor','none')
addParameter(p,'meancolor','k')
addParameter(p,'plotminmax',true)
addParameter(p,'minmaxcolor',[.7 .7 .7])
addParameter(p,'minmaxedgecolor','none')
addParameter(p,'minmaxfacealpha',1)

parse(p,SW,varargin{:});

% parameters
distadjust = p.Results.distadjust;
left       = p.Results.left;
right      = p.Results.right;
boundedl   = p.Results.boundedline;

if ~left && ~right
    warning('Empty SWATHobj: nothing to plot')
    return;
end


if ~left
    ny = ceil(length(SW.disty)/2)+1;
    SW.Z = SW.Z(ny:end,:);
elseif ~right
    ny = floor(length(SW.disty)/2);
    SW.Z = SW.Z(1:ny,:);
end

z_min  = min(SW.Z,[],1);
z_max  = max(SW.Z,[],1);
z_mean = mean(SW.Z,1,'omitnan')';
z_std  = std(SW.Z,0,1,'omitnan')';
dist   = SW.distx+distadjust;

% draw SWATHobj
ax = gca;
if ishold(ax)
    ISHOLD = true;
else
    ISHOLD = false;
    cla(ax, 'reset')
    hold(ax,'on')
    box(ax,'on')
end


if (boundedl) %exist('boundedline','file') && 
    if p.Results.plotminmax
    hp(3) = patch([dist(:); flipud(dist(:))],...
                  [z_min(:); flipud(z_max(:))],...
                  p.Results.minmaxfacecolor,...
                  'EdgeColor',p.Results.minmaxedgecolor,...
                  'FaceAlpha',p.Results.minmaxfacealpha);
    end
    hp(2) = patch([dist(:); flipud(dist(:))],...
                  [z_mean(:)-z_std(:); flipud(z_mean(:)+z_std(:))],...
                  p.Results.facecolor,...
                  'EdgeColor',p.Results.edgecolor,...
                  'FaceAlpha',p.Results.facealpha);
    

    hp(1) = plot(dist,z_mean,'color',p.Results.meancolor);

else
    hp(1) = plot(dist,z_mean,'-','Color',p.Results.meancolor); 
    hp(2) = plot([dist;nan;dist],[z_mean-z_std;nan;z_mean+z_std],'-','color',clr);
    if p.Results.plotminmax
        hp(3) = plot([dist;nan;dist],[z_min nan z_max],':',...
            'color',clr);
    end
end

drawnow
xlabel(sprintf('Distance along profile (%s)',SW.xyunit))
ylabel(sprintf('Z (%s)',SW.zunit))
if p.Results.plotminmax
    legend(hp,{'Mean','+/- St.Dev.','Min/Max'})
else
    legend(hp,{'Mean','+/- St.Dev.'})
end

if ~ISHOLD
    hold(ax,'off')
end

if nargout == 1
    h = hp;
end

