function DEM = elevateminima(DEM,maxarea)

%ELEVATEMINIMA elevate regional minima in a DEM to their lowest neighbor
%
% Syntax
%
%     DEMe = elevateminima(DEM)
%     DEMe = elevateminima(DEM,area)
%
% Description
%
%     Digital elevation models often feature spurious pits, e.g.
%     erroneous depressions that usually comprise only one to a few pixel
%     with lower elevations than their neighbors. The function fillsinks
%     is usually used to deal with topographic sinks, but does not
%     account for nested pits and thus cannot handle spurious pits located
%     inside larger true sinks or sinks that should be carved rather than
%     filled. elevateminima identifies regionalminima using the IPT
%     function imregionalmin and elevates them to their highest surrounding
%     neighbor. A maximum area (default = 1 pixel) can be defined to
%     control which local maximima to elevate.
%
%     Note that regional minima "are connected components of pixels with a 
%     *constant* intensity value, and whose external boundary pixels all 
%     have a higher value" (see documention of imregionalmin).
%
% Input arguments
%
%     DEM       digital elevation model (GRIDobj)
%     maxarea   maximum area of regional minima in pixels to be elevated
%               (default = 1)
%
% Output arguments
%
%     DEMe      digital elevation model (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEM2 = elevateminima(DEM,2);
%     imageschs(DEM,DEM2-DEM);
%
%
% See also: GRIDobj/fillsinks, imregionalmin, imreconstruct
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 30. January, 2013

narginchk(1,2);
if nargin == 1;
    maxarea = 1;
else
    validateattributes(maxarea,{'numeric'},{'scalar','integer','>=',1},'elevateminima','maxarea',2);
end



INAN = isnan(DEM.Z); 
DEM.Z(INAN) = inf;
I = imregionalmin(DEM.Z,8);
I(INAN) = true;
I = imclearborder(I);
I = xor(bwareaopen(I,maxarea+1),I);

IX = find(I);
nrrows = DEM.size(1);
IXoffset = [-1 -1+nrrows nrrows 1+nrrows 1 1-nrrows -nrrows -1-nrrows];

DEM.Z(I) = inf;

for n = 1:8;
    DEM.Z(IX) = min(DEM.Z(IX),DEM.Z(IX+IXoffset(n)));
end

if maxarea > 1
    STATS = regionprops(I,DEM.Z,'MinIntensity','PixelIdxList');
    for r = 1:numel(STATS);
        DEM.Z(STATS(r).PixelIdxList) = STATS(r).MinIntensity;
    end
end

DEM.Z(INAN) = nan;
    
