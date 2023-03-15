function [x,y,varargout] = STREAMobj2XY(S,varargin)

%STREAMOBJ2XY convert instance of STREAMobj to NaN-separated X and Y coordinates
%
% Syntax
%
%     [x,y] = STREAMobj2XY(S)
%     [x,y,a,b,...] = STREAMobj2XY(S,A,B,...)
%
% Description
%
%     STREAMobj2XY returns the NaN-punctuated vectors x and y which can be
%     used for easy plotting using the plot function. With additional input
%     arguments that are instances of GRIDobj or node attribute lists,
%     additional vectors are generated with the respective values of the
%     grids A, B, etc. at the locations x and y.
%
% Input arguments
%
%     S       streams (class STREAMobj)
%     A,B,... grids (class GRIDobj) or node attributes (e.g. as returned by
%             the function STREAMobj/streamorder
%
% Output arguments
%
%     x,y     coordinate vectors
%     a,b,... grid values at locations specified by x and y.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     [x,y,d] = STREAMobj2XY(S,S.distance);
%     plot(x,y)
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 30. January, 2013


order = S.orderednanlist;
I     = ~isnan(order);
x     = nan(size(order));
x(I)  = S.x(order(I));
y     = nan(size(order));
y(I)  = S.y(order(I));

nrnodes = numel(S.x);
if nargin > 1
    varargout = cell(numel(varargin),1);
    for r = 1:numel(varargin)
        
        if isa(varargin{r},'GRIDobj')
            % extract values from GRIDobj
            validatealignment(S,varargin{r})        
            varargout{r}    = nan(size(order));
            varargout{r}(I) = double(varargin{r}.Z(S.IXgrid(order(I))));
        elseif numel(varargin{r}) == nrnodes
            % extract values from node attribute list
            varargout{r}    = nan(size(order));
            varargout{r}(I) = double(varargin{r}(order(I)));
        end

    end
    
end
    
    
    