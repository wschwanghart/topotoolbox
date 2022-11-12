function R = meltonruggedness(FD,DEM)
%MELTONRUGGEDNESS Melton ruggedness
%
% Syntax
%
%     R = meltonruggedness(FD,DEM)
%
% Description
%    
%     Melton ruggedness number (MNR) is a simple flow accumulation related
%     index, calculated as difference between maximum and minimum elevation
%     in catchment area divided by square root of catchment area size. The
%     calculation is performed for each grid cell, therefore minimum
%     elevation is same as elevation at cell's position. Due to the
%     discrete character of a single maximum elevation, flow calculation is
%     simply done with the D8 flow direction.
%     (Description from SAGA GIS, see link in source)
%     
%
% Input arguments
%
%     FD   FLOWobj
%     DEM  Digital elevation model
%
% Output arguments
%
%     R    Melton ruggedness
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     R = meltonruggedness(FD,DEM);
%     imagesc(R)
%     imageschs(DEM,R)
%
% See also: GRIDobj/roughness, FLOWobj/upslopestats, FLOWobj/flowacc
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 12. November, 2022


% link
% https://sourceforge.net/p/saga-gis/code/ci/master/tree/saga-gis/src/tools/terrain_analysis/ta_hydrology/melton_ruggedness.cpp

A = flowacc(FD)*DEM.cellsize^2;
M = upslopestats(FD,DEM,'max');

R = (M-DEM)/sqrt(A);


