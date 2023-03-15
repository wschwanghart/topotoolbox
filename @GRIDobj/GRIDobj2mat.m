function [A,x,y] = GRIDobj2mat(DEM)

%GRIDobj2mat convert GRIDobj to matrix and coordinate vectors
%
% Syntax
%
%     [dem,X,Y] = GRIDobj2mat(DEM)
%
% Description
%
%     convert GRIDobj to matrix and coordinate vectors 
%
% Input
%
%     DEM       instance of GRIDobj class
% 
% Output
%
%     dem       matrix 
%     x,y       coordinate vectors
%                  
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     [Z,x,y] = GRIDobj2mat(DEM);
%     plot(x,Z(20,:))
%     xlabel('x-coordinate [m]')
%     ylabel('Elevation [m]')
%
% See also: GRIDobj 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017

A = DEM.Z;
[x,y] = refmat2XY(DEM.refmat,DEM.size);
