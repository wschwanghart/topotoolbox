function extractconncomps(S)

%EXTRACTCONNCOMPS interactive stream network selection
%
% Syntax
%
%     extractconncomps(S)
%
% Description
%
%     extractconncomps displays a figure and plots the stream network S.
%     Here you can mouse-select individual connected components and export
%     them to the workspace.
%
% Input arguments
%
%     S    instance of a STREAMobj
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     S  = STREAMobj(FD,A>1000);
%     extractconncomps(S)
%
% See also: STREAMobj/klargestconncomps
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017 


CS = STREAMobj2cell(S);
fx = figure;

% Create push button
btn1 = uicontrol('Style', 'pushbutton', 'String', 'Export to workspace',...
        'Position', [20 20 150 20],...
        'Callback', @exporttoworkspace);    
btn1 = uicontrol('Style', 'pushbutton', 'String', 'Switch selection',...
        'Position', [170 20 150 20],...
        'Callback', @switchselection);    


ax = axes('parent',fx);
hold on
nc = numel(CS);
clr_nonsel = [.6 .7 .7];
clr_sel    = [0 0 0];
sel = false(nc,1);
for r = 1:nc
    h(r) = plot(CS{r},'color',clr_nonsel);
    h(r).ButtonDownFcn = @highlight;
end
box on
axis image
title('Click network to select/unselect.')


%% --- Highlight function ----
function highlight(h,~)

if isequal(h.Color,clr_nonsel)
    h.Color = clr_sel;
else
    h.Color = clr_nonsel;
end
end

%% --- Switch selection ----
function switchselection(hh,~)
    for r2 = 1:numel(h)
        if isequal(h(r2).Color,clr_sel)
            h(r2).Color = clr_nonsel;
        else
            h(r2).Color = clr_sel;
        end
    end
end

%% --- Export to workspace ----
function exporttoworkspace(hh,~)

    I = false(nc,1);
    for r2 = 1:numel(h);
        I(r2) = isequal(h(r2).Color,clr_sel);
    end
    n = nnz(I);
    if n==0
        warndlg('No streams available for export.');
        return
    elseif n == 1
        S2 = CS{I};
    else
        CSS = CS(I);
        S2 = union(CSS{:});
    end

    prompt = {'Enter variable name:'};
    ptitle = 'Export';
    plines = 1;
    pdef = {'S'};
    answer = inputdlg(prompt, ptitle, plines, pdef);
    if ~isempty(answer) && isvarname(answer{1})
        assignin('base',answer{1},S2);
    else
        return
    end


    end

end




   
% 
% 
% nrc = numel(S.x);
% M = sparse(double(S.ix),double(S.ixc),true,nrc,nrc);
% 
% [L,nc] = conncomps(S);
% 
% 
% if nargin == 1;
% 
% hFig = figure('Units','normalized','OuterPosition',[0 0 1 1]);
% hAx  = axes('Parent',hFig);
% cmap = jet(nc);
% 
% for r = 1:nc;
%     M2 = spdiags(L==r,0,nrc,nrc)*double(M);
%     [x,y] = gplot(M2,[S.x S.y]);
%     outlet = find((sum(M2,1)' > 0) & (sum(M2,2) == 0));
%     
%     plot(hAx,x,y,'Color',cmap(r,:));
%     hold on
%     plot(hAx,S.x(outlet),S.y(outlet),'*','Color',cmap(r,:));
%     text(S.x(outlet),S.y(outlet),[' ' num2str(r)]);
% end
% 
% axis image
% set(hAx,'Xtick',[],'Ytick',[]);
% 
% prompt = {'Enter comma separated list of component indices to extract:'};
% 
% dlg_title = 'Extract connected components';
% num_lines = 1;
% options.Resize='on';
% answer = inputdlg(prompt,dlg_title,num_lines,{''},options);
% 
% % evaluate answer
% c = regexp(answer,'(?:\d*)','match');
% c = str2double(c{1});
% else
%     c = cc;
%     if c > nc;
%         error('TopoToolbox:wronginput','The supplied number exceeds the number of components');
%     end
% end
% 
% % adapt new STREAMobj to the reduced network
% L     = ismember(L,c);
% I     = L(S.ix);
% S.ix  = S.ix(I);
% S.ixc = S.ixc(I);
% 
% IX    = cumsum(L);
% S.ix  = IX(S.ix);
% S.ixc = IX(S.ixc);
% 
% S.x   = S.x(L);
% S.y   = S.y(L);
% S.IXgrid   = S.IXgrid(L);
