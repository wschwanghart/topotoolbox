function extout = setextent(extent,ax)

%SETEXTENT set current axis extent 
%
% Syntax
%
%     setextent(extent)
%     setextent(DEM)
%     setextent(S)
%     setextent(P)
%     setextent(xy)
%	  setextent(MS)
%     setextent(...,ax)
%     ext = setextent(...)
%     
% Description
%
%     getextent and setextent are two small wrapper functions (around set
%     and get) to quickly zoom to a specified extent in an axis object.
%
% Input parameters
%
%     extent    cell array as returned by getextent
%     DEM       GRIDobj
%     S         STREAMobj
%     P         PPS
%     xy        nx2 array with n coordinate pairs
%     MS        mapstruct (must have the fields X and Y)
%     ax        target axes (default is gca)
%
% Example
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,'minarea',1000);
%     S   = modify(S,'distance',20000);
%     S   = klargestconncomps(S);
%     imageschs(DEM);
%     hold on
%     plot(S,'w')
%     hold off
%     setextent(S)
%
% See also: IMAGESCHS, GETEXTENT, PADEXTENT
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 23. January, 2020

if nargin == 1
    ax = gca;
end

zoom reset
if iscell(extent)
    set(ax,{'xlim','ylim'},extent)
elseif isa(extent,'GRIDobj')
    [x,y] = getoutline(extent);
    set(ax,'xlim',[min(x) max(x)],'ylim',[min(y) max(y)]);
elseif isstruct(extent)
    x = [extent.X];
    y = [extent.Y];
    set(ax,'xlim',[min(x) max(x)],'ylim',[min(y) max(y)]);
elseif isa(extent,'STREAMobj')
    v = info(extent,'boundingbox');
    set(ax,'xlim',v(1:2),'ylim',v(3:4));
elseif isa(extent,'PPS')
	v = info(extent.S,'boundingbox');
    set(ax,'xlim',v(1:2),'ylim',v(3:4));
elseif isnumeric(extent)
    x = extent(:,1);
    y = extent(:,2);
	
    set(ax,'xlim',[min(x) max(x)],'ylim',[min(y) max(y)]);
    
end

if nargout == 1
    extout = getextent(ax);
end