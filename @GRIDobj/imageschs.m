function rgb = imageschs(DEM,A,varargin)

%IMAGESCHS plot hillshade image with overlay
%
% Syntax
%
%     imageschs(DEM)
%     imageschs(DEM,A)
%     imageschs(DEM,A,pn,pv,...)
%     RGB = imageschs(...)
%
% Description
%
%     Hillshading is a very powerful tool for relief depiction. imageschs
%     calculates a hillshade of a digital elevation model and colors it
%     according to a second grid or matrix. DEM must be an instance of  
%     GRIDobj and A must be a GRIDobj or matrix. 
%
%     imageschs allows to set a number of parameter name/value pairs that
%     control lighting direction and representation of missing values
%     (nan).
%
%     If called with an output variable, imageschs does not plot the
%     hillshade but returns an RGB image. The hillshading algorithm follows
%     the logarithmic approach to shaded relief representation of Katzil
%     and Doytsher (2003).
%
%     Some of the parameter settings (e.g. colorbarylabel) require MATLAB's
%     graphics engine introduced in R2014b and will throw errors when used
%     with older versions.
%
% Input
%
%     DEM         digital elevation model (GRIDobj)
%     A           coloring matrix or GRIDobj
%
% Parameter name/value pairs
%
%     caxis            two element vector defining the value range. Default 
%                      is [min(A) max(A)].  
%     colorbar         false or true (default)
%     colorbarlabel    string. Title for the colorbar
%     colorbarylabel   string. Label for the y-axis of the colorbar
%     colormap         string for colormap name or [ncol x 3] matrix. Note 
%                      that if NaNs or Infs are found in A, the colormap  
%                      must not have more than 255 colors. Default: 'jet'
%     percentclip      scalar prc (%) that truncates the displayed range to  
%                      the prc's and 100%-prc's percentile of the data in A.
%                      This parameter is ignored if 'caxis' is defined.
%                      Default is prc=0.
%     truecolor        three element vector (rgb) with values between 0 and 1  
%                      that indicates how true values are plotted if A is 
%                      logical.
%                      Default is [0 1 0].
%     falsecolor       three element vector (rgb) with values between 0 and 1  
%                      that indicates how false values are plotted if A is 
%                      logical.
%                      Default is [1 1 1].
%     nancolor         three element vector (rgb) with values between 0 and 1  
%                      that indicates how NaNs and Infs are plotted 
%                      Default is [1 1 1].
%     brighten         Scalar between -1 and 1. Shift intensities of all 
%                      colors in the colormap. If <0, then colormap will be
%                      darkened, and if >0, it will be brightened.
%                      Default is 0
%     usepermanent     controls whether the hillshade is retained in 
%                      memory as persistent variable. Default is false. If 
%                      set to true, and if the DEM has the same size as the
%                      persistently stored hillshade, the function will
%                      reuse this variable, thus avoiding to recalculate
%                      hillshading.
%     medfilt          use median filter to smooth hillshading 
%                      (default=false)
%     method           'surfnorm' or 'mdow'. mdow is the multidirectional
%                      oblique hillshade algorithm.
%     azimuth          azimuth angle of illumination, (default=315)
%     altitude         altitude angle of illumination, (default=60)
%     exaggerate       elevation exaggeration (default=2). Increase to
%                      pronounce elevation differences in flat terrain. If
%                      the DEM is in geographic coordinates, use a value of
%                      0.000003.
%     ticklabels       'default', 'nice' or 'none'
%     tickstokm        true or {false}. If set to true, coordinates will be
%                      divided by 1000.
%     gridmarkers      two element vector with [dx dy] spacing of + markers
%     gridmarkercolor  three element vector (rgb) or color abbreviations
%                      as given in LineSpec (default = 'k')
%                 
%
% Output
%
%     RGB         [DEM.size 3] image (UINT8) with values between 0 and 255
%
%
% Example 1: Hillshade of DEM with colors derived from elevations. 
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     imageschs(DEM);
%
% Example 2: Hillshade of DEM with non-default colormap 
%     
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     imageschs(DEM,[],'colormap',landcolor);
%
% Example 3: Gray-colored hillshading only
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     imageschs(DEM,[],'colormap',[.9 .9 .9],'colorbar',false);
%
% Example 4: Slope map with underlying hillshade
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     G   = gradient8(DEM);
%     imageschs(DEM,G);
%
% Example 5: Nice axis tick labels
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     imageschs(DEM,DEM,'ticklabels','nice','colorbar',true);
% 
% Example 6: Map logical GRIDobj with custom colors
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     L   = DEM>1500;
%     imageschs(DEM,L,'ticklabels','nice','colorbar',true,...
%                       'truecolor',[1 0 0],'falsecolor',[0 0 1]);
%
% Example 7: Calculate hillshade RGB and export as geotiff
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     RGB = imageschs(DEM);
%     geotiffwrite('file.tif',RGB,DEM.georef.SpatialRef,...
%            'GeoKeyDirectoryTag',DEM.georef.GeoKeyDirectoryTag);
%
% Example 8: Calculate hillshade RGB and display using imshow
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     RGB = imageschs(DEM);
%     RGB = flipud(RGB);
%     [~,R] = GRIDobj2im(DEM);
%     imshow(RGB,R)
%     axis xy
%
% References
%
%     Katzil, Y., Doytsher, Y. (2003): A logarithmic and sub-pixel approach
%     to shaded relief representation. Computers & Geosciences, 29,
%     1137-1142.
%
% See also: HILLSHADE
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 13. June, 2016


% Change log
%
% 23.1.2014: line 198: A = gray2ind(mat2gray(A,ncolors); 
%            replaced with
%            A = gray2ind(mat2gray(A,alims),ncolors);
% 02.6.2015: changed help and updated to R2014b
% 04.3.2016: added option medfilt and added example
% 18.3.2016: added option percentclip
% 13.6.2016: added option makepermanent
% 14.6.2016: added options 
% 06.4.2018: new example

persistent H

narginchk(1,inf);
nargoutchk(0,1);

% if A is not supplied to the function, coloring will be according to
% values in DEM
if nargin == 1 || (nargin>=2 && isempty(A))
    A = DEM;
end

%% Default values can be changed here
% -----------------------------------
defaultcolormap = 'parula';
defaulttruecolor = [0 1 0]; 
defaultfalsecolor = [1 1 1]; 
defaultnancolor = [1 1 1];
defaultexaggerate = 1;
defaultazimuth  = 315;
defaultaltitude = 60;
defaultcolorbar = true;
defaultmedfilt  = false;
% -----------------------------------
%%

% Parse inputs
p = inputParser;
p.FunctionName = 'GRIDobj/imageschs';
% required
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'A',@(x) isa(x,'GRIDobj') || ismatrix(A));
% optional
addParamValue(p,'colormap',defaultcolormap,@(x)(ischar(x) || size(x,2)==3));
addParamValue(p,'caxis',[],@(x) numel(x) == 2);
addParamValue(p,'percentclip',[],@(x) isscalar(x) && x>=0 && x<50);
addParamValue(p,'truecolor',defaulttruecolor,@(x) isequal(size(x),[1 3]) && (max(x)<=1) && (min(x) >= 0));
addParamValue(p,'falsecolor',defaultfalsecolor,@(x) isequal(size(x),[1 3]) && (max(x(:))<=1) && (min(x(:)) >= 0));
addParamValue(p,'nancolor',defaultnancolor,@(x) isequal(size(x),[1 3]) && (max(x(:))<=1) && (min(x(:)) >= 0));
addParamValue(p,'exaggerate',defaultexaggerate,@(x) isscalar(x) && x>0);
addParamValue(p,'azimuth',defaultazimuth ,@(x) isscalar(x) && x>0);
addParamValue(p,'altitude',defaultaltitude,@(x) isscalar(x) && x>0);
addParamValue(p,'colorbar',defaultcolorbar,@(x) isscalar(x));
addParamValue(p,'medfilt',defaultmedfilt,@(x) isscalar(x));
addParamValue(p,'ticklabels','default',@(x) ischar(x));
addParamValue(p,'gridmarkers',[],@(x) numel(x) == 1 || numel(x) == 2);
addParamValue(p,'gridmarkercolor','k');
addParamValue(p,'useparallel',true);
addParamValue(p,'usepermanent',false);
addParamValue(p,'colorbarlabel',[],@(x) ischar(x));
addParamValue(p,'colorbarylabel',[],@(x) ischar(x));
addParamValue(p,'tickstokm',false,@(x) isscalar(x));
addParamValue(p,'method','surfnorm');
addParamValue(p,'brighten',0,@(x) x>=-1 & x <= 1);
parse(p,DEM,A,varargin{:});

% required
DEM        = p.Results.DEM;
A          = p.Results.A;
colmapfun  = p.Results.colormap;
cbar       = p.Results.colorbar;
truecol    = p.Results.truecolor;
falsecol   = p.Results.falsecolor;
exag       = p.Results.exaggerate;
azi        = p.Results.azimuth;
alti       = p.Results.altitude;
nancolor   = p.Results.nancolor;
usepermanent  = p.Results.usepermanent;
tokm           = p.Results.tickstokm > 0;
colorBarLabel  =p.Results.colorbarlabel;
colorBarYLabel =p.Results.colorbarylabel;
meth       = validatestring(p.Results.method,{'default','surfnorm','mdow'});


ticklabels = validatestring(p.Results.ticklabels,{'default','none','nice'});
gridmarkers= p.Results.gridmarkers;
gridmarkercolor = p.Results.gridmarkercolor;

% check if input matrices align
validatealignment(DEM,A)

if isa(A,'GRIDobj')
    A = A.Z;
end

% constrain color range to values given in caxis
if ~isempty(p.Results.caxis)
    A(A<p.Results.caxis(1)) = p.Results.caxis(1);
    A(A>p.Results.caxis(2)) = p.Results.caxis(2);
elseif ~isempty(p.Results.percentclip)
    qclip = p.Results.percentclip/100;
    [n,edges] = histcounts(A(~isnan(A(:))),'Normalization','cdf');
    lval = edges(find(n>=qclip,1,'first'));
    uval = edges(find(n<(1-qclip),1,'last'));
    if lval == uval
        warning('TopoToolbox:imageschs','percent clip returns flat matrix');
        A(:,:) = lval;
    else
        A = max(A,lval);
        A = min(A,uval);
    end
        
end

% coordinate matrices
[x,y] = refmat2XY(DEM.refmat,DEM.size);

% convert coordinates to km if wanted
if tokm    
    x = x*1e-3;
    y = y*1e-3;
end

% nr of colors
nhs = 256;

% calculate hillshading
if usepermanent && isequal(size(H),DEM.size)

else
    H = hillshade(DEM,'exaggerate',exag,'azimuth',azi,'altitude',alti,'useparallel',p.Results.useparallel,'method',meth);
%     switch meth
%         case 'default'    
%             H = hillshade(DEM,'exaggerate',exag,'azimuth',azi,'altitude',alti,'useparallel',p.Results.useparallel);
%         case 'mdow'
%             H = hillshademdow(DEM,'exaggerate',exag,'useparallel',p.Results.useparallel);
%     end
    H = H.Z;
    Inan = isnan(H);
    if any(Inan(:))
        H(Inan) = 1;
        clear Inan
    else
        clear Inan
    end
    
    % median filtering, if required
    if p.Results.medfilt
        H = medfilt2(H,[3 3],'symmetric');
    end
    
    H = gray2ind(H,nhs);
end

% derive coloring
if ~isa(A,'logical')
    Inan = isnan(A(:)) | isinf(A(:));
    if any(Inan(:))
        nans = true;
        A(Inan) = nan;
    else
        nans = false;
        clear Inan
    end
    
    if isa(colmapfun,'char')        
        ncolors = 256-nans;
        colmapfun = str2func(lower(colmapfun));
        cmap = colmapfun(ncolors);
    else
        ncolors = size(colmapfun,1);
        if nans && (ncolors >= 256)
            error('TopoToolbox:GRIDobj',['NaNs found in the second input grid. \n'...
                  'Please provide colormap with less than 256 colors']);
        else
            cmap = colmapfun;
        end        
    end
    
    if cbar && isempty(p.Results.caxis)
        alims = [min(A(:)) max(A(:))];
    elseif cbar && ~isempty(p.Results.caxis)
        alims = sort(p.Results.caxis,'ascend');
    else
        alims = [min(A(:)) max(A(:))];
    end
    
    A = gray2ind(mat2gray(A,double(alims)),ncolors);
    
else
    ncolors = 2;
    cmap = [falsecol; truecol];
	nans = false;
    alims = [0 1];
end
    
% create colormap for indexing
cmap = cmap(:);
cmap = bsxfun(@times,cmap,linspace(0,1,nhs));
cmap = reshape(cmap,[ncolors 3 nhs]);
cmap = permute(cmap,[3 1 2]);
cmap = reshape(cmap,[ncolors*nhs 3]);

if p.Results.brighten
    cmap = brighten(cmap,p.Results.brighten);
end

% create image that indexes into the new colormap
IND  = uint16(H+1) + nhs*uint16(A) + 1;

% handle NaNs
if nans
    cmapnan   = bsxfun(@times,nancolor,linspace(0,1,nhs)');
    IND(Inan) = uint16(H(Inan)) + nhs*(ncolors) +1;% unclear if this is ok...
    cmap      = [cmap;cmapnan];
end

% same as ind2rgb but returns a mxnx3 matrix with uint8 data
cmapUINT8 = uint8(round(cmap*256));
% see Rob Campbell's mat2im
% http://www.mathworks.de/matlabcentral/fileexchange/26322-mat2im
RGB=reshape(cmapUINT8(IND(:),:),[size(IND),3]);

% use permanent
if ~usepermanent
    H = [];
end

% plot
if nargout == 0
    imagesc(x,y,RGB);
    axis xy
    axis image
    
    % add colorbar if needed
    if cbar
        if alims(1) ~= alims(2) 
            caxis(alims);
        end
        colormap(gca,cmap(nhs:nhs:nhs*ncolors,:));
        cc = colorbar;%('location','south');
        if ~isempty(colorBarLabel)
            title(cc,colorBarLabel);
        end
        if ~isempty(colorBarYLabel)
            ylabel(cc,colorBarYLabel);
        end
    end
    
    % plot nice ticklabels if wanted
    switch ticklabels
        case 'none'
            set(gca,'XTickLabel',{},'YTickLabel',{});
        case 'nice'
            xticklocs = get(gca,'XTick');
            yticklocs = get(gca,'YTick');
            
            set(gca,'XTick',xticklocs([1 end]))
            set(gca,'YTick',yticklocs([1 end]))
            
            set(gca,'XTickLabel',num2str(xticklocs([1 end])','%d'));
            set(gca,'YTickLabel',num2str(yticklocs([1 end])','%d'));
            
            % rotate tick labels if matlab 2014b or newer available
            if ~verLessThan('matlab','8.4')
                set(gca,'YTickLabelRotation',90);
            end
                       
    end
    
    % plot grid
    if ~isempty(gridmarkers)
        if numel(gridmarkers) == 1
            gridmarkers = [gridmarkers gridmarkers];
        end
          
        xgridmarkers = unique(x-rem(x,gridmarkers(1)));
        ygridmarkers = unique(y-rem(y,gridmarkers(2)));
        hold on
        [xx,yy] = meshgrid(xgridmarkers,ygridmarkers);
        plot(xx(:),yy(:),'+','Color',gridmarkercolor);
        hold off
    end
        
    
elseif nargout == 1
    rgb = RGB;
end
