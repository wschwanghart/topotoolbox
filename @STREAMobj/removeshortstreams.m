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

% downstream distance
dsdistance = zeros(numel(S.x),1);
for r=1:numel(S.ix);
    dsdistance(S.ixc(r)) = max(dsdistance(S.ixc(r)),...
        dsdistance(S.ix(r)) + ...
        sqrt((S.x(S.ixc(r))-S.x(S.ix(r)))^2 + (S.y(S.ixc(r))-S.y(S.ix(r)))^2));
end

keepvertices = true(numel(S.x),1);
keepvertices(dsdistance<=d & s==1) = false;

for r=numel(S.ix):-1:1;
    if (s(S.ixc(r)) == 1) && (s(S.ix(r)) == 1);
        keepvertices(S.ix(r)) = keepvertices(S.ixc(r));
    end
end

% adapt new STREAMobj to the reduced network
L     = keepvertices;
I     = L(S.ix);
S.ix  = S.ix(I);
S.ixc = S.ixc(I);

IX    = cumsum(L);
S.ix  = IX(S.ix);
S.ixc = IX(S.ixc);

S.x   = S.x(L);
S.y   = S.y(L);
S.IXgrid   = S.IXgrid(L);

S = clean(S);










