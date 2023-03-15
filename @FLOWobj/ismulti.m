function tf = ismulti(FD,typetest)
%ISMULTI check if FD is multi or single flow direction
%
% Syntax
%
%     tf = ismulti(FD)
%
% Description
%
%     ISMULTI returns true if a FLOWobj contains multiple flow directions
%
% Input arguments
%
%     FD     FLOWobj
%
% Output arguments
%
%     tf     true or false
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     ismulti(FD)
% 
%     ans =
% 
%       logical
% 
%        0
%
% See also: FLOWobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. March, 2016

if nargin == 1
    typetest = true;
end

if typetest
tf = strcmp(FD.type,'multi');
else 
tf = any(histcounts(FD.ix,1:(prod(FD.size)+1))>1);
end
