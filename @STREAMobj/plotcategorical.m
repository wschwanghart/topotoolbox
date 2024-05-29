function [h,hcb] = plotcategorical(S,c,varargin)

%PLOTCATEGORICAL Plot categorical variable along a single river
%
% Syntax
%
%     h = plotcategorical(S,c)
%
% Description
%
%     This function plots a categorical variable c along a single river
%     stored in the STREAMobj S. c can be a node-attribute list or a
%     GRIDobj. By default, the function will plot patches of contiguous
%     values along the upstream distance of the river. Y-values are set to
%     0 and 1, by default. 
%     
% Input arguments
%
%     S     STREAMobj (single stream only)
%     c     GRIDobj or node-attribute values with categorical values. c can
%           be numeric, logical, or categorical. 
%
%     Parameter name/value pairs
%
%     'colormap'    nc*3 array with RGB values. nc is the number of unique 
%                   values in c. Alternatively, you can provide an
%                   anonymous function (e.g. @parula) or a name of a 
%                   function (e.g. 'parula') that creates a colormap.
%     'colorbar'    true or {false}. If colorbar is true, then a colorbar
%                   will be plotted. Note that axes colormap will then be
%                   adjusted and the color range is set to 0.5 to nc+0.5,
%                   where nc is the number of unique values in c.
%     'ztop'        scalar or node-attribute list. Default = 1. The values
%                   determine the upper Y-values of the patches.
%     'zbot'        scalar or node-attribute list. Default = 0. The values
%                   determine the lower Y-values of the patches.
%     'distance'    node-attribute list with distance values. By default, 
%                   the distance is S.distance.
%     'dunit'       distance unit. The default is 'm'. If set to 'km',
%                   distance values are divided by 1000.
%     'EdgeColor'   default = 'none'
%     'FaceAlpha'   default = 1
%     'alphagrad'   false. If true, then the patches will transition from
%                   fully opaque at the top to fully transparent at the
%                   bottom.
%     'parent'      default = gca
%     
% Output arguments
%
%     h     handles of patches
%     hcb   handle to colorbar
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S  = STREAMobj(FD, minarea = 1000);
%     S  = trunk(klargestconncomps(S));
%     c  = mod(S.distance,10000)>5000;
%     plotcategorical(S,c,'zbot',0,'ztop',getnal(S,DEM))
%     hold on
%     plotdz(S,DEM);
%
% See also: STREAMobj/plotdz, STREAMobj/plotc
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 24. May, 2014

p = inputParser;
p.FunctionName = 'STREAMobj/plotcategorical';

addRequired(p,'S',@(x) isa(x,'STREAMobj') && info(S,'nrchan')==1)
addRequired(p,'c',@(x) isa(x,'GRIDobj') || isnal(S,x))
addParameter(p,'colormap',@parula)
addParameter(p,'colorbar',false)
addParameter(p,'distance',S.distance)
addParameter(p,'dunit','m')
addParameter(p,'EdgeColor','none')
addParameter(p,'FaceAlpha',1)
addParameter(p,'parent',gca)
addParameter(p,'alphagrad',false)
addParameter(p,'ztop',1,@(x) isnal(S,x) || isscalar(x))
addParameter(p,'zbot',0,@(x) isnal(S,x) || isscalar(x))

parse(p,S,c,varargin{:})

d = p.Results.distance;
ax = p.Results.parent;
if ~ishold(ax)
    ax = newplot(ax);
end

switch lower(p.Results.dunit)
    case 'm'
        dunit = p.Results.dunit;
    case 'km'
        dunit = p.Results.dunit;
        d = d/1000;
    otherwise
        warning('Distance unit not known.')
        dunit = '-';
end

if isnal(S,c)
    % do nothing
else
    c = getnal(S,c);
end

if ~isnal(S,p.Results.ztop)
    ztop = getnal(S)+p.Results.ztop;
    zvar = false;
else
    ztop = p.Results.ztop;
    if numel(unique(ztop)) == 0
        zvar = false;
    else
        zvar = true;
    end
end
if ~isnal(S,p.Results.zbot)
    zbot = getnal(S)+p.Results.zbot;
    zvar2 = false;
else
    zbot = p.Results.zbot;
    if numel(unique(zbot)) == 0
        zvar2 = false;
    else
        zvar2 = true;
    end
end
zvar = zvar || zvar2;

if iscategorical(c)
    catnames = categories(c);
    [~,~,c] = unique(c);
else
    [catnames,~,c] = unique(c);
    catnames = cellfun(@num2str,num2cell(catnames),'UniformOutput',false);
end
nClasses = numel(catnames);

[~,~,d,c,zbot,ztop] = STREAMobj2XY(S,d,c,zbot,ztop);
d = d(1:end-1); % remove trailing nans
c = c(1:end-1);
ztop = ztop(1:end-1);
zbot = zbot(1:end-1);

% Get categories and their colors


% Function 
if isa(p.Results.colormap,"function_handle")
    clr = p.Results.colormap(nClasses);
elseif ischar(p.Results.colormap)
    clr = str2func(p.Results.colormap);
    clr = clr(nClasses);
else
    clr = p.Results.colormap;
    validateattributes(clr,{'numeric'},{'nrows',nClasses,'ncols',3,'>=',0,'<=',1})
end

% draw patches for each contiguous region of classes
I = [true; diff(c) ~= 0];
I(end) = true;
ix = find(I);

% FaceVertexAlphaData
if ~(p.Results.alphagrad)
    fa  = p.Results.FaceAlpha;
    fag = [];
    usefag = false;
else
    fa = 'interp';
    fag = [];
    usefag = true;
end


if ~zvar
    for r = 1:numel(ix)-1
        ixl = ix(r);
        ixr = ix(r+1);
        X   = d([ixl ixr ixr ixl ixl]);
        Y   = [zbot([ixl ixr]); ztop([ixr ixl]); zbot(ixl)];
        if usefag
            fag = [0 0 1 1 0]';
        end
        hout(r) = patch(ax,X(:),Y(:),clr(c(ixl),:),...
            'FaceAlpha',fa,...
            'EdgeColor',p.Results.EdgeColor,...
            'FaceVertexAlphaData',fag);
    end
else
    for r = 1:numel(ix)-1
        ixl = ix(r);
        ixr = ix(r+1);
        ixx = ixl:ixr;
        X   = d([ixx ixx(end:-1:1) ixx(1)]);
        Y   = [zbot(ixx); ztop(ixx(end:-1:1)); zbot(ixx(1))];
        if usefag
            fag = [zeros(numel(ixx),1);ones(numel(ixx),1);0];
        
        end
        hout(r) = patch(ax,X(:),Y(:),clr(c(ixl),:),...
            'FaceAlpha',fa,...
            'EdgeColor',p.Results.EdgeColor,...
            'FaceVertexAlphaData',fag);
        if usefag
            set(hout(r),'FaceAlpha','interp')
        end
    end
end

% The axes may lie behind the patches. The next line brings the axes to the
% front.
set(ax,'Layer','Top')

% Make colorbar
if p.Results.colorbar
    colormap(ax,clr)
    clim(ax,[0.5 max(c)+0.5]);
    hcbout = colorbar(ax,"eastoutside",'Xtick',1:max(c),...
        'Xticklabel',catnames);
else
    hcbout = [];
end


if nargout >= 1
    h = hout;
end
if nargout == 2
    hcb = hcbout;
end


