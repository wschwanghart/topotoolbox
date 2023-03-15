function DEM = filter(DEM,varargin)

%FILTER 2D-filtering of DEMs with different kernels 
%
% Syntax
%
%     DEMF = filter(DEM)
%     DEMF = filter(DEM,method)
%     DEMF = filter(DEM,method,kernelsize)
%
% Description
%
%     The function filter is a wrapper around various image filtering
%     algorithms including mean, sobel, median etc. So far, only filters
%     with rectangular kernels are supported.
%
% Input
%
%     DEM         digital elevation model (GRIDobj)
%     method      'mean' (default), 'sobel', 'scharr', 'median', 'wiener'
%     kernelsize  size of moving window, default [3 3]
%
% Output
%
%     DEMF        filtered digital elevation model (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEMf = filter(DEM,'wiener',[11 11]);
%     imageschs(DEM,DEM-DEMf)
%     
% 
% See also: CONV2, FILTER2, MEDFILT2, WIENER2
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 26. August, 2020



validmethods = {'mean','average','median',...
                'sobel','scharr','wiener',...
                'std'};

% Parse inputs
p = inputParser;
p.StructExpand  = true;
p.KeepUnmatched = false;
p.FunctionName = 'filter'; 
addRequired(p,'DEM',@(x) isequal(class(x),'GRIDobj'));
addOptional(p,'method','mean',@(x) ischar(x));
addOptional(p,'kernel',[3 3],@(x) isscalar(x) || (numel(x) == 2) || isempty(x));
parse(p,DEM,varargin{:});

DEM     = p.Results.DEM;
method  = validatestring(p.Results.method,validmethods);
ws      = p.Results.kernel;
if isscalar(ws)
    ws = [ws ws];
end
dem    = DEM.Z;

% pad DEM if there are NaNs
inan    = isnan(dem);
flagnan = any(inan(:));
if flagnan
    [~,L]   = bwdist(~inan); 
    dem     = dem(L);
end

switch method
    case {'mean','average'}
        padsize = ceil(ws/2);
        dem     = padarray(dem,padsize,'replicate');
        W = ones(ws);
        dem = conv2(dem,W./sum(W(:)),'same');
        dem = dem(padsize(1)+1:end-padsize(1),padsize(2)+1:end-padsize(2));
    case 'median'
        dem     = medfilt2(dem,ws,'symmetric');
    case {'sobel','scharr'}
        if any(ws~=3)  
            warning('TopoToolbox:GRIDobj',...
                ['The method ' method ' only works with a 3x3 kernel']);
        end
        switch method
            case 'sobel'
                ky = [1 2 1; 0 0 0; -1 -2 -1];
            case 'scharr'
                ky = [3 10 3; 0 0 0; -3 -10 -3];
        end
        
        kx = ky';
        padsize = ceil(ws/2);
        dem     = padarray(dem,padsize,'replicate');
        dem     = hypot(conv2(dem,ky,'valid'),conv2(dem,kx,'valid'));
        dem     = dem(2:end-1,2:end-1);
        
    case 'wiener'
        dem     = wiener2(dem,ws);
    case 'std'
        dem     = stdfilt(dem,true(p.Results.kernel(1),p.Results.kernel(2)));
end

% and set nans at their previous position
if flagnan
    dem(inan) = nan;
end

% write output
DEM.Z = dem;
DEM.name = [method ' filter'];



