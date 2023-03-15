function ht = surf(DEM,varargin)

%SURF surface plot for GRIDobj
%
% Syntax
%
%     surf(DEM)
%     surf(DEM,A)
%     surf(...,pn,pv,...)
%     h = ...
%
% Description
%
%     surf for GRIDobj overloads the surf command and thus provides fast
%     access to 3D visualization of digital elevation models. Note that
%     GRIDobj/surf automatically resamples the DEM so that the maximum of
%     rows or columns does not exceed 1000. This ensures that the surface
%     is efficiently drawn.
%
% Input arguments
%
%     DEM          digital elevation model (GRIDobj)
%     A            grid to define color (GRIDobj) 
%
% Parameter name/values
%
%     'exaggerate'   height exaggeration, default = 1
%      
%     and all property name/value pairs allowed by surface objects.
%
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     surf(DEM)
%     camlight
%
% Example 2
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     surf((DEM-800)*2,'block',true,'baselevel',-900*2,'sea',true,'exaggerate',1); camlight; colormap(ttcmap((DEM-800)*2,'cmap','france'))
%     ax = gca;
%     ax.Clipping = 'off';
%     axis off
%     ax.Projection = 'perspective';
%
%     FD = FLOWobj(DEM);
%     S  = STREAMobj(FD,'minarea',1000);
%     S  = modify(S,'upstreamto',(DEM-800)>0);
%     hold on
%     plot3(S,DEM-800,'b');
%
% See also: imageschs
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 30. January, 2013



%% Parse inputs
p = inputParser;
p.KeepUnmatched = true;
addOptional(p,'A',[]);
addParameter(p,'maxrowsorcols',2000)
addParameter(p,'exaggerate',1);
addParameter(p,'block',false);
addParameter(p,'baselevel',[]);
addParameter(p,'sea',false);
addParameter(p,'sealevel',0);
addParameter(p,'seaalpha',0.5);

% Parse
parse(p,varargin{:});

% create variables
maxrowsorcols = p.Results.maxrowsorcols;
exagg = p.Results.exaggerate;
block = p.Results.block;
baselevel = p.Results.baselevel;

if block
    minz = min(DEM);
    maxz = max(DEM);
    if isempty(baselevel)
    baselevel = minz-(maxz-minz)*0.2;
    end
end    
sea   = p.Results.sea;
sealevel = p.Results.sealevel;
seaalpha = p.Results.seaalpha;
%%

if max(DEM.size)>maxrowsorcols
    
    resamplecellsize = max(DEM.size)/(maxrowsorcols-1) * DEM.cellsize;
    DEM = resample(DEM,resamplecellsize);
    
    if ~isempty(p.Results.A)
        A = p.Results.A;
        A = resample(A,DEM);
        overlay = true;
    else
        overlay = false;
    end
         
    
else
    if ~isempty(p.Results.A)
        A = p.Results.A;
        validatealignment(DEM,A);
        overlay = true;
    else
        overlay = false;
    end
end


[x,y] = refmat2XY(DEM.refmat,DEM.size);

% Create cellarray from p.Unmatched
pn = fieldnames(p.Unmatched);
pv = struct2cell(p.Unmatched);

pnpv = [pn pv];
pnpv = pnpv';
pnpv = pnpv(:)';

if overlay
    h = surf(x,y,double(DEM.Z),double(A.Z),pnpv{:});
else
    h = surf(x,y,double(DEM.Z),pnpv{:});
end

exaggerate(gca,exagg);
shading interp
% camlight

if block 
    
    if any(isnan(DEM))
        error('TopoToolbox:wronginput','DEM must not have NaNs')
    end
    facecolor  = [.6 .6 .6];
    facecolordark = brighten(facecolor,-0.2);
    
    xp = x(:);
    xp = [xp(1); xp; xp(end:-1:1)];
    yp = repmat(y(1),size(xp));
    zp = DEM.Z(1,:);
    zp = zp(:);
    zp = [baselevel;zp;repmat(baselevel,size(zp))];
    pa(1) = patch(xp,yp,zp,facecolor,'EdgeColor','none','FaceColor',facecolor);
    
    yp = repmat(y(end),size(xp));
    zp = DEM.Z(end,:);
    zp = zp(:);
    zp = [baselevel;zp;repmat(baselevel,size(zp))];
    pa(2) = patch(xp,yp,zp,facecolor,'EdgeColor','none','FaceColor',facecolor);
    
    yp = y(:);
    yp = [yp(1); yp; yp(end:-1:1)];
    xp = repmat(x(1),size(yp));
    zp = DEM.Z(:,1);
    zp = zp(:);
    zp = [baselevel;zp;repmat(baselevel,size(zp))];
    pa(3) = patch(xp,yp,zp,facecolordark,'EdgeColor','none','FaceColor',facecolordark);
    
    xp = repmat(x(end),size(yp));
    zp = DEM.Z(:,end);
    zp = zp(:);
    zp = [baselevel;zp;repmat(baselevel,size(zp))];
    pa(4) = patch(xp,yp,zp,facecolordark,'EdgeColor','none','FaceColor',facecolordark);
    
    set(pa,'FaceLighting','none')
    
    
    axislim = axis;
    if baselevel < min(DEM)
        axislim(5) = baselevel;
    end
    axis(axislim)
end
        

if sea 
    
    if any(isnan(DEM))
        error('TopoToolbox:wronginput','DEM must not have NaNs')
    end
    facecolor  = [1 1 1];
    facecolordark = [0.8 .8 .8];
    
    xp = x(:);
    xp = [xp(1); xp; xp(end:-1:1)];
    yp = repmat(y(1),size(xp));
    zp = DEM.Z(1,:);
    zp = min(zp(:),0);
    zp = double([sealevel;zp;repmat(sealevel,size(zp))]);
    pa(1) = patch(xp,yp,zp,facecolor,'EdgeColor','none','FaceColor',facecolor,'FaceAlpha',seaalpha);
    
    yp = repmat(y(end),size(xp));
    zp = DEM.Z(end,:);
    zp = min(zp(:),0);
    zp = double([sealevel;zp;repmat(sealevel,size(zp))]);
    pa(2) = patch(xp,yp,zp,facecolor,'EdgeColor','none','FaceColor',facecolor,'FaceAlpha',seaalpha);
    
    yp = y(:);
    yp = [yp(1); yp; yp(end:-1:1)];
    xp = repmat(x(1),size(yp));
    zp = DEM.Z(:,1);
    zp = min(zp(:),0);
    zp = double([sealevel;zp;repmat(sealevel,size(zp))]);
    pa(3) = patch(xp,yp,zp,facecolordark,'EdgeColor','none','FaceColor',facecolordark,'FaceAlpha',seaalpha);
    
    xp = repmat(x(end),size(yp));
    zp = DEM.Z(:,end);
    zp = min(zp(:),0);
    zp = double([sealevel;zp;repmat(sealevel,size(zp))]);
    pa(4) = patch(xp,yp,zp,facecolordark,'EdgeColor','none','FaceColor',facecolordark,'FaceAlpha',seaalpha);
    
    I  = DEM.Z>0;
%     H  = zeros(DEM.size);
%     H(I) = nan;
%     % Texturemap
%     I  = ~I;
%     TM = repmat(I,1,1,3);

    hold on
    h = surf(x,y,+I,'FaceAlpha',0.5,'FaceColor',facecolor,'EdgeColor','none');
    hold off
    
    
    set(pa,'FaceLighting','none')
    
    if ~isempty(baselevel)
    axislim = axis;
    axislim(5) = baselevel;
    axis(axislim)
    end
    
    
    

end
    
if nargout == 1
    ht = h;
end


function exaggerate(axes_handle,exagfactor)

% elevation exaggeration in a 3D surface plot
%
% Syntax
%
%     exaggerate(axes_handle,exagfactor)
%
% Description
%
%     exaggerate is a simple wrapper for calling set(gca...). It controls
%     the data aspect ratio in a 3D plot and enables elevation
%     exaggeration.   
%
% Input
%
%     axes_handle   digital elevation model
%     exagfactor    exaggeration factor (default = 1)
%
% Example
%
%     load exampledem
%     for r = 1:4;
%     subplot(2,2,r);
%     surf(X,Y,dem); exaggerate(gca,r);
%     title(['exaggeration factor = ' num2str(r)]);
%     end
%
% 
% See also: SURF
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 6. September, 2010

if nargin == 1
    exagfactor = 1;
end
axis(axes_handle,'image');
set(axes_handle,'DataAspectRatio',[1 1 1/exagfactor]);



