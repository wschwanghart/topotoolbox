function s = streamorder(S,type)

%STREAMORDER calculate Strahler Stream Order from STREAMobj
%
% Syntax
%
%     s = streamorder(S)
%     streamorder(S)
%
% Description
%
%     The Strahler Stream Order is a way to classify rivers based on a
%     hierarchy of tributaries. First-order streams don't have tributaries.
%     When two first-order streams come together, they form a second-order
%     stream. When two second-order streams come together, they form a
%     third-order stream. When two streams of different order confluence,
%     they form a stream of maximum order of both. streamorder returns the
%     Strahler Stream Order based on an instance of the STREAMobj (S). 
%
%     s = streamorder(S)     returns a column vector with node attributes 
%     of the directed graph represented by STREAMobj.
%
%
% Input
%   
%     S         STREAMobj     
%     type      {'Strahler'} or 'shreve'
%
% Output
%
%     s         stream order for each node in STREAMobj
%
% Example
%     
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     so = streamorder(S);
%     plotc(S,so)
%     colorbar
% 
% See also: STREAMobj, FLOWobj/streamorder, STREAMobj/plotstreamorder
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. January, 2013


if nargin == 1;
    type = 'strahler';
else
    type = validatestring(type,{'strahler','shreve'},'STREAMobj/streamorder','type',2);
end

switch type
    case 'strahler'
        nrc = numel(S.x);
        s   = zeros(nrc,1);
        
        for r = 1:numel(S.ix);
            if s(S.ix(r)) == 0;
                s(S.ix(r)) = 1;
            end
            
            if s(S.ix(r)) < s(S.ixc(r))
            elseif s(S.ix(r)) == s(S.ixc(r))
                s(S.ixc(r)) = s(S.ixc(r))+1;
            else
                s(S.ixc(r)) = s(S.ix(r));
            end
        end
        
    case 'shreve'
        s = streampoi(S,'channelheads','logical');
        s = double(s);
        
        for r = 1:numel(S.ix);
            s(S.ixc(r)) = s(S.ix(r))+s(S.ixc(r));
        end
end

  
    
    
