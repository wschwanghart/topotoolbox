function cmap = ttscm(name,n,percrange)

%TTSCM Scientific colormaps by Fabio Crameri
%
% Syntax
%     
%     cmap = ttscm(name)
%     cmap = ttscm(name,n)
%     cmap = ttscm(name,percrange)
%     ttscm
%     allowedcmaps = ttscm
%
% Description
%
%     TTSCM provides access to scientific colormaps compiled by Fabio 
%     Crameri (http://www.fabiocrameri.ch/colourmaps.php). 
%
%     TTSCM(name) returns a 255*3 matrix with rgb values.
%
%     TTSCM(name,n) returns a n*3 matrix with rgb values.
%
%     TTSCM(name,n,percrange) enables to define a percentile range of the
%     total colormap to be used. A scalar > 50 thereby defines the central
%     percentile range, whereas a scalar < 50 defines the upper and lower
%     limits of the percentile range. percrange can be a two element vector
%     with the upper and lower percentile bounds of the colormap range to
%     be used.
%     
%     TTSCM without input and output arguments shows a figure with
%     available colormaps.
%
%     TTSCM without input arguments and one output argument returns the
%     names of the allowed colormaps as cell array.
%
% Input arguments
%
%     name           name of colormap (run ttscm without in- and output
%                    arguments for list)
%     n              number of colors in colormap    
%     percrange      percentile range of the colormap to be used (see
%                    description above)
%     
% Output arguments
%
%     cmap           n*3 colormap
%     allowedcmaps   allowed colormaps
%
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     C   = curvature(DEM);
%     lims = prcclip(C,2,true);
%     clr = ttscm('vik');
%     imageschs(DEM,C,'colormap',clr,'caxis',lims);
%
% Example 2: Show colormaps using the Big Tujunga Catchment data
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     cmaps = ttscm;
%     tiledlayout('flow','TileSpacing','none');
%     for r = 1:numel(cmaps); 
%        nexttile 
%        imageschs(DEM,[],'colormap', ttscm(cmaps{r}),...
%                         'colorbar',false,'ticklabels','none',...
%                         'useperm',true);
%        title(cmaps{r});
%     end
%
%     
% References: 
%
%      Crameri, F., (2018). Scientific colour-maps. Zenodo. 
%      http://doi.org/10.5281/zenodo.4491293
%
%      Crameri, F. (2018), Geodynamic diagnostics, scientific visualisation 
%      and StagLab 3.0, Geosci. Model Dev., 11, 2541-2562. 
%      doi:10.5194/gmd-11-2541-2018
%
%
% See also: GRIDobj/IMAGESCHS, ttcmap
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 10. June, 2021

allowedcmaps = {'acton' 'bamako', 'batlow', 'hawaii', 'imola' 'nuuk' ...
                'buda' ...
                'devon','davos','oslo','bilbao','lajolla',...
                'grayC','broc','cork','vik','lisbon','tofino',...
				'berlin','turku','tokyo','lapaz','roma','oleron', ...
                'brocO','corkO','vikO', 'bam', 'bamO', 'batlowK', 'batlowW',...
				'bukavu','fes'};

				
% get location of this function
p = fileparts(mfilename('fullpath'));
if nargin == 0 && nargout == 0
    figure('MenuBar','none','Color','w','Toolbar','none',...
           'Name', 'Scientific Colour Maps', 'NumberTitle', 'off');
           
	imshow([p filesep 'private' ...
            filesep '+ScientificColourMaps_FabioCrameri.png'],...
            'InitialMagnification','fit')
	return
elseif nargin == 0 && nargout == 1
    cmap = sort(allowedcmaps);
    return
end

cmaptype = validatestring(name,allowedcmaps);
if nargin == 1
    n = 255;
else
    validateattributes(n,{'numeric'},{'>',1},'ttscm','n',2)
end

cmap = loadcmap(cmaptype);

if nargin == 3
    if isscalar(percrange)
        if percrange > 100 || percrange < 0
            error('TopoToolbox:ttscm','Wrong percentile range')
        end
        
        if percrange >= 50
            percrange = [50-percrange/2 50+percrange/2];
        else
            percrange = [percrange 100-percrange];
        end
        
    elseif numel(percrange) == 2
        
        percrange = sort(percrange,'ascend');
        if any(percrange > 100) || any(percrange < 0)
            error('TopoToolbox:ttscm','Wrong percentile range')
        end
        
    end
    ncolors = size(cmap,1);
    cmap = cmap((round(percrange(1)/100*ncolors)+1):round(percrange(2)/100*ncolors),:);
end

ncolors = size(cmap,1);
if n ~= ncolors
	cmap    = interp1((1:ncolors)',cmap,linspace(1,ncolors,n));
end



