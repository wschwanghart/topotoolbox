function varargout = as(P,outtype)

%AS Convert PPS object into various data formats 
%
% Syntax
%
%     S = as(P,'STREAMobj')
%     G = as(P,'GRIDobj')
%     [MSS,MSP] = as(P,'mapstruct')
%     [MSS,MSP] = as(P,'mapshape')
%     [xs,ys,xp,yp] = as(P,'nanpunctuated')
%     [G,p,xdata,ydata] = as(P,'graph') or as(P,'digraph')
%
% Description
% 
%     AS converts a PPS object P into a different format that can be
%     worked with in MATLAB or exported to other software.
%
%     S = as(P,'STREAMobj') extracts the STREAMobj that underlies the PPS
%     object. This is the same as S = P.S
%
%     G = as(P,'GRIDobj') generates a GRIDobj. G has values of 1 along
%     streams, values of 2 at points, and zeros elsewhere.
%
%     [MSS,MSP] = as(P,'mapstruct') creates structure arrays that can be
%     exported to shapefiles (see also function PPS/shapewrite). MSS is a
%     structure array containing the stream network and MSP contains the
%     points.
%
%     [MSS,MSP] = as(P,'mapshape'). Same as 'mapstruct', but the resulting
%     variables are mapshape objects.
%
%     [xs,ys,xp,yp] = as(P,'nanpunctuated') creates nan-punctuated vectors
%     of the stream network xs and ys, and vectors of the point coordinates
%     xp and yp. The data can be used to easily customize plotting of a PPS
%     object.
%
%     [G,p,xdata,ydata] = as(P,'graph') or as(P,'digraph') creates a graph
%     and digraph from the stream network. p is the node id of the points
%     on the network. xdata and ydata are vectors with the coordinates of
%     the node vertices.
%
% Input arguments
%
%     P      PPS/object
%     
% Output arguments
%
%     see description above
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,'minarea',10000);
%     P   = PPS(S,'rpois',0.0005);
%     [G,p,xdata,ydata] = as(P,'digraph');
%     plot(G,'-','xdata',xdata,'ydata',ydata,...
%         'arrowsize',0,'marker','none')
% 
% See also: GRIDobj, FLOWobj, STREAMobj, PPS/shapewrite,
%           STREAMobj/STREAMobj2mapstruct, STREAMobj/STREAMobj2XY
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 10. January, 2020

switch lower(outtype)
    case 'streamobj'
        varargout{1} = P.S;
    case 'gridobj'
        G = GRIDobj(P.S);
        G.Z(P.S.IXgrid) = 1;
        G.Z(P.S.IXgrid(P.PP)) = 2;
        varargout{1} = G;
    case 'mapstruct'
        varargout{1} = STREAMobj2mapstruct(P.S);
        
        if nargout == 2
            [x,y] = points(P);
            id = (1:numel(x))';
            varargout{2} =  struct('Geometry','Point',...
            'X',num2cell(x),...
            'Y',num2cell(y),...
            'id',num2cell(id));
            
        end
    case 'mapshape'
        [MSS,MSP] = as(P,'mapstruct');
        varargout{1} = mapshape(MSS);
        varargout{2} = mapshape(MSP);
    case 'nanpunctuated'
        [varargout{1},varargout{2}] = STREAMobj2XY(P.S);
        [varargout{3},varargout{4}] = points(P);
    case {'graph' 'digraph'}
        xdata   = P.S.x;
        ydata   = P.S.y;
        weights = hypot(xdata(P.S.ix)-xdata(P.S.ixc),...
                        ydata(P.S.ix)-ydata(P.S.ixc));
        M = sparse(P.S.ix,P.S.ixc,weights,numel(xdata),numel(ydata)); 
        T = table((1:numel(P.S.x))',points(P,'nal'),'VariableNames',{'pts','ispoint'});
        
        switch lower(outtype)
            case 'graph'
                G = graph(M+M',T);
            case 'digraph'
                G = digraph(M,T);
        end
        varargout{1} = G;
        varargout{2} = P.PP;
        varargout{3} = xdata;
        varargout{4} = ydata;
    
end