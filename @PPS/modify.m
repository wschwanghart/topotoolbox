function [P,nalix,markix] = modify(P,varargin)

%MODIFY modify instance of PPS to meet user-defined criteria
%
% Syntax
%
%      P2 = modify(P,'pn',pv)
%      [P2,nalix,markix] = ...
%
% Description
%
%      PPS/modify modifies the underlying stream network of the PPS object
%      P. The function is a wrapper for the STREAMobj/modify function but
%      differs in that it carries the points stored in P along. For a
%      detailed description, see the help of STREAMobj/modify.
%
%      help STREAMobj/modify
%
% Input arguments
%
%      P        Instance of PPS
%      'pn',pv  Parameter value pair (see STREAMobj/modify)
%
% Output arguments
%
%      P2       new instance of PPS
%      nalix    linear index into node-attribute list underlying P
%      markix   linear index into marks associated with the point pattern
%               in P
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     S = removeshortstreams(S,100);
%     S = clean(S);
%     P = PPS(S,'rpois',0.001,'z',DEM);
%
%     subplot(1,2,1)
%     plot(P)
%     
%     P2 = modify(P,'streamorder','>2');
%     subplot(1,2,2)
%     plot(P2)
%
% See also: STREAMobj/modify
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. May, 2021

% Apply STREAMobj/modify
[S,nalix] = modify(P.S,varargin{:});

% Elevations are stored as node attribute list
if ~isempty(P.z)
    z = P.z(nalix);
else
    z = [];
end

% Index into markers
I = ismember(P.S.IXgrid(P.PP),P.S.IXgrid(nalix));
markix = find(I);

% Finally, create new PPS instance
P = PPS(S,'PP',P.S.IXgrid(P.PP(I)));
% Add elevation
P.z = z;