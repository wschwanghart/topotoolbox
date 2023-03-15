function h = plot(S,varargin)

%PLOT plot instance of STREAMobj
%
% Syntax
%
%     plot(S)
%     plot(S,...)
%     h = ...;
%
% Description
%
%     plot overloads the builtin plot command and graphically outputs the
%     stream network in S. Additional plot options can be set in the same
%     way as with the builtin plot command.
%
% Input arguments
%
%     S
%
% Parameter name/value pairs
%
%     'cdc'   {false} or true. If true, a customized data cursor will be 
%             constructed so that upstream distance will be displayed in 
%             addition to x and y coordinates. However, this will result in
%             an error when data cursor is used to underlying images.
%    
%     'labeldist' place distance indicators on the stream network.
%             Value can be either a numeric array with distance values or a
%             cell array. The cell array must have two elements where the 
%             first element contains a numeric array with distances and the
%             second element contains a string with a line specification
%             syntax (e.g., '+r' will draw ret crosses). If you specify
%             'text' as second element, than text labels will be drawn at
%             the distance locations. When called with one output argument,
%             h will contain several handles to line and text objects.
%
%     other pn/pv pairs same as plot command
%
% Output arguments
%
%     h       handle to line object
%
% Example
%     
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     S  = STREAMobj(FD,'minarea',1000);
%     plot(S)
% 
%     % plot with distance indicators at 10000 and 50000 m from outlets
%     % with default distance indicators
%     plot(S,'labeldist',[10000 50000])
%     % with red squares
%     plot(S,'labeldist',{[10000 50000] 'sr'})
%     % or with text labels
%     plot(S,'labeldist',{[10000 30000] 'text'})
%
%
% See also: STREAMobj, STREAMobj/plotdz
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 30. January, 2013
        

% go through varargin to find 'exaggerate'
TF = strcmpi('cdc',varargin);
ix = find(TF);
if ~isempty(ix)
    cdc = varargin{ix+1};
    varargin([ix ix+1]) = [];
else
    cdc = false;
end

TF = strcmpi('labeldist',varargin);
ix = find(TF);
if ~isempty(ix)
    labeldist = varargin{ix+1};
    varargin([ix ix+1]) = [];
else
    labeldist = [];
end


if cdc || ~isempty(labeldist)
    [x,y,d] = STREAMobj2XY(S,S.distance);
else
    [x,y] = STREAMobj2XY(S);
end

ht = plot(x,y,varargin{:});

%%%
if cdc
    dcm_obj = datacursormode();
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,d})
end
%%%
if ~isempty(labeldist)
    xyd = [];
    % find nan separators
    ixnansep = [0; find(isnan(d))];
    if iscell(labeldist);
        linestyle = labeldist{2};
        labeldist = labeldist{1};
    else
        linestyle = '+b';
    end
    labeldist = double(labeldist(:));
    for r = 1:numel(ixnansep)-1;
        ix = ixnansep(r)+1:ixnansep(r+1)-1;
        xy = interp1(d(ix),[x(ix) y(ix)],labeldist,'linear',nan);
        inan = ~isnan(xy(:,1));
        if any(inan);
            xyd = [xyd;[xy(inan,:) labeldist(inan)]];            
        end
    end
    ih = ishold;
    hold on;
    if ~strcmpi(linestyle,'text');
        ht(2) = plot(xyd(:,1),xyd(:,2),linestyle);
    else
        h(2) = plot(xyd(:,1),xyd(:,2),'.k');
        htext = text(xyd(:,1),xyd(:,2),num2str(xyd(:,3)),'verticalalignment','baseline');
        h    = [h htext(:)'];
    end
    if ~ih;
        hold off;
    end
end

if nargout == 1
    h = ht;
end



function txt = myupdatefcn(~,event_obj,d)
% Customizes text of data tips
pos = get(event_obj,'Position');
try
I   = get(event_obj, 'DataIndex');
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['Distance: ',num2str(d(I))]};
catch
end
    

