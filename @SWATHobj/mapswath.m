function [OUT] = mapswath(SW,M,varargin)
%MAPSWATH obtains Z values along a SWATHobj from an arbitrary GRIDobj
%
% Syntax
%
%     [OUT] = mapswath(SW,M)
%     [OUT] = mapswath(SW,M,method)
%
% Description
%
%     MAPSWATH interpolates Z values from the GRIDobj M, at the data point
%     locations of the SWATHobj SW. If SW and M are of different spatial
%     reference (projection), MAPSWATH tries to use the mapping toolbox to 
%     convert between units.
%
% Input arguments
%
%     SW      instance of SWATHobj
%     M       instance of GRIDobj
%     method  interpolation method, can be any of 'linear','nearest',
%             'spline','pchip','cubic'; by default 'linear'
%
% Output arguments
%
%     OUT    instance of SWATHobj
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     SW = SWATHobj(DEM,'dx',200,'dy',200);
%     G = gradient8(DEM,'degree');
%     SWG = mapswath(SW,G);
%     plotdz(SWG)
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.edu)
% Date: May, 2015


narginchk(2,3)
if nargin == 2;
    method = 'linear';
else
    method = validatestring(varargin{1},...
        {'linear','nearest','spline','pchip','cubic'},'mapswath','method',3);
end

% check inputs
if ~isa(SW,'SWATHobj')
    error('First input hast to be of class SWATHobj');
    
elseif ~isa(M,'GRIDobj')
    error('Second input hast to be of class GRIDobj');
    
else
    OUT = SW;
    OUT.Z(:) = nan;
    ix = find(~isnan(SW.Z));
    OUT.Z(ix) = interp(M,OUT.X(ix),OUT.Y(ix),method);
    
    % check results
    if any(isnan(OUT.Z))
        warning('Output SWATHobj contains NaN. Input arguments may not fully overlap.')
    elseif all(isnan(OUT.Z))
        warning('Output SWATHobj contains only NaN. Input arguments may not overlap.')
    end

    OUT.name = {'SWATHobj created from SWATHobj:';SW.name};
    OUT.georef = M.georef;
    
end








