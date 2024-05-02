function S = removeshortstreams(S,d)

%REMOVESHORTSTREAMS Remove first order streams with a length less than specified
%
% Syntax
% 
%     S2 = removeshortstreams(S,d)
%
% Description
%
%     Digital stream networks sometimes include short first order streams
%     that should not be included in further analysis. removeshortstreams
%     enables to remove such dangling first order streams the length in  
%     map units of which is less or equal to d.
%
% Input arguments
%
%     S     streams (class STREAMobj)
%     d     length of first order streams in map units to be removed.
%
% Output arguments
%
%     S2    streams (class STREAMobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S2 = removeshortstreams(S,2000);
%     plot(S)
%     hold on
%     plot(S2)     
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017

narginchk(2,2)

validateattributes(d,{'numeric'},{'scalar','>',0},'removeshortstreams','d',2);

% calculate streamorder
s = streamorder(S);
I = s == 1;

S2 = subgraph(S,I);
c  = aggregate(S2,S2.distance,'method','drainagebasins','aggfun',@max);
c  = c>d;
I  = nal2nal(S,S2,c,true);
S  = subgraph(S,I);

