classdef subplotlabel < handle
    
%SUBPLOTLABEL Label subplots (works only for 2D plots so far)
%
% Syntax
%
%     h = subplotlabel(ax,letter)
%     h = subplotlabel(fig,letter)
%     h = subplotlabel(...,pn,pv,...)
%
% Description
%
%     subplotlabel adds labels to the corner of each panel of a composite
%     figure.
%
%     subplotlabel(ax,'c') assigns the letter 'c' to the axis ax.
%
%     subplotlabel(fig,'A') assigns the letters 'A','B',... to the axes in
%     the figure with the handle fig. letter must be 'A','a', 'I', 'i', 
%     or '1'.
%     
% Input arguments
%     
%     ax     axes handle (e.g. gca)
%     fig    figure handle (e.g. gcf)
%
%     Parameter name/value pairs
% 
%     'Location'    one of the following values: {'northwest'}, 'southwest', 
%                   'northeast', 'southeast', 'northwestoutside'
%     'FontSize'    14
%     'FontWeight'  {'normal'} or 'bold'
%     'FontAngle'   {'normal'} or 'italic'
%     'Color'       Font color. {'k'}
%     'BackgroundColor'   {'none'} or any other way to define colors
%
%     Applicable if called with fig handle only
%
%     'order'       direction of labelling {'rightdown'} or 'downright'
%
% Output arguments
%
%     h    instance of subplotlabel
%
% Methods for class subplotlabel
%
%     bigger(h,s)       increases the font size by s. Default is 1.
%     smaller(h,s)      decreases the font size by s. Default is 1.
%     resetposition(h)  resets the position. 
%     upper(h)          sets all labels to upper case
%     lower(h)          sets all labels to lower case
%     italic(h)         sets all labels to italic
%     bold(h)           sets all labels to bold
%     normal(h)         sets all labels to normal
%     delete(h)         removes all labels
%
% Example 1
%
%     for r = 1:6; subplot(2,3,r); plot(rand(5)); end
%     h = subplotlabel(gcf,'a','location','southwest','order','downright');
%
% Example 2 (requires TopoToolbox)
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     ax = subplot(1,2,1);
%     imageschs(DEM)
%     h = subplotlabel(ax,'A','location','northwest','color','w');
%     ax(2) = subplot(1,2,2);
%     imageschs(DEM,gradient8(DEM))
%     h(2) = subplotlabel(ax(2),'B','location','southeast','color','w');
%     bigger(h,5)
%
% See also: text
%
% Note: roman numericals are computed using Francois Beauducel function
% NUM2ROMAN (https://github.com/beaudu/romanum)
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 22. February, 2021

    properties
        label
        LIMIX
    end
    
    methods
        function h = subplotlabel(varargin)
            
            
            p = inputParser;
            addRequired(p,'ax')
            addRequired(p,'value')
            addParameter(p,'location','northwest');
            addParameter(p,'order','rightdown');
            addParameter(p,'FontSize',14);
            addParameter(p,'Margin',1);
            addParameter(p,'EdgeColor','none');
            addParameter(p,'BackgroundColor','none');
            addParameter(p,'FontWeight','normal');
            addParameter(p,'FontAngle','normal');
            addParameter(p,'Color','k');
            % addParameter(p,'Clipping','on');
            parse(p,varargin{:})
            
            if isa(p.Results.ax,'matlab.ui.Figure')
                hax = findall(gcf,'Type','axes');
                
                if isempty(hax)
                    error('Figure must contain at least one axis')
                end
                
                switch p.Results.value
                    case 'A'
                        letters = ('A':'Z')';
                        letters = cellstr(letters);
                    case 'a'
                        letters = ('a':'z')';
                        letters = cellstr(letters);
                    case '1'
                        letters = num2cell(1:numel(hax));
                        letters = cellfun(@num2str,letters,'UniformOutput',false);
                    case {'i' 'I'}
                        letters = cellfun(@roman,num2cell(1:numel(hax)),'UniformOutput',false);
                        switch p.Results.value
                            case 'i'
                                letters = lower(letters);
                        end
                        
                    otherwise
                        error('TopoToolbox:subplotlabel',...
                            ['If the first argument is a handle to a figure, letter \newline' ...
                             'must be ''A'', ''a'', ''I'' ''i'' or ''1''.']);
                end
                
                % Sort axes from top-left to bottom-right
                pos = zeros(numel(hax),2);
                % Left upper corner of axes
                for r = 1:numel(hax)
                    pos(r,1) = hax(r).Position(1);
                    pos(r,2) = sum(hax(r).Position([2 4])); 
                end
 
                switch lower(p.Results.order)
                    case 'rightdown'
                        [~,ix] = sortrows(pos,[-2 1]);
                    case 'downright'
                        [~,ix] = sortrows(pos,[1 -2]);
                    otherwise 
                        error('unknown ordering')
                end   
               
                hax = hax(ix);
                Results = p.Results;
                Results = rmfield(Results,{'ax','value'});
                
                for r = 1:numel(hax)
                    hh(r) = subplotlabel(hax(r),letters{r},Results);
                end
                
                h = hh;
                return
            else
                ax = p.Results.ax;
                letter = p.Results.value;
                
                xl = [0 1]; %xlim(ax);
                yl = [0 1]; %ylim(ax);
                switch lower(p.Results.location)
                    case {'northeast','ne'}
                        IXX = 2;
                        IXY = 2;
                        valign = 'top';
                        halign = 'right';
                    case {'northwest','nw'}
                        IXX = 1;
                        IXY = 2;
                        valign = 'top';
                        halign = 'left';
                    case {'southwest','sw'}
                        IXX = 1;
                        IXY = 1;
                        valign = 'bottom';
                        halign = 'left';
                    case {'southeast','se'}
                        IXX = 2;
                        IXY = 1;
                        valign = 'bottom';
                        halign = 'right';
                    case {'northwestoutside'}
                        IXX = 1;
                        IXY = 2;
                        valign = 'bottom';
                        halign = 'left';
                end
                
                locx = xl(IXX);
                locy = yl(IXY);
                
                switch halign
                    case 'right'
                        letter = [' ' letter ' '];
                    case 'left'
                        letter = [' ' letter ' '];
                end
                
                if ishold(ax)
                    ihold = true;
                else
                    ihold = false;
                end
                
                hold(ax,'on');
                h.label = text(ax,locx,locy,letter,...
                    'VerticalAlignment',valign,...
                    'HorizontalAlignment',halign,...
                    'Color',p.Results.Color,...
                    'FontSize',p.Results.FontSize,...
                    'Margin',p.Results.Margin,...
                    'BackgroundColor',p.Results.BackgroundColor,...
                    'EdgeColor',p.Results.EdgeColor,...
                    'FontAngle',p.Results.FontAngle,...
                    'FontWeight',p.Results.FontWeight,...
                    'Units','normalized',...
                    'Clipping','on');
                
                
                if ~ihold
                    hold(ax,'off');
                end
                
                h.LIMIX  = [IXX, IXY];
            
            end
        end
        function resetposition(h)
            % Set subplot labels to initial position
            for r = 1:numel(h)
                xl = [0 1]; 
                yl = [0 1];
                h(r).label.Position = [xl(h(r).LIMIX(1)) yl(h(r).LIMIX(2))];
            end
        end
        function upper(h)
            % Upper case letters
            for r = 1:numel(h)
                set(h(r).label,'String',upper(h(r).label.String));
            end
        end
        function lower(h)
            % Lower case letters
            for r = 1:numel(h)
                set(h(r).label,'String',lower(h(r).label.String));
            end
        end
        function bigger(h,val)
            % Increase font size
            if nargin == 1
                val = 1;
            end
            for r = 1:numel(h)
                h(r).label.FontSize = h(r).label.FontSize + val;
            end
        end
        function smaller(h,val)
            % Decrease font size
            if nargin == 1
                val = 1;
            end
            for r = 1:numel(h)
                h(r).label.FontSize = h(r).label.FontSize - val;
            end
        end
        function italic(h)
            % Italic font
            for r = 1:numel(h)
                h(r).label.FontAngle = 'italic';
            end
        end
        function bold(h)
            % Bold font
            for r = 1:numel(h)
                h(r).label.FontWeight = 'bold';
            end
        end
        function normal(h)
            % Normal font
            for r = 1:numel(h)
                h(r).label.FontWeight = 'normal';
                h(r).label.FontAngle  = 'normal';
            end
        end
        function delete(h)
            % Delete subplot labels
            
            for r = 1:numel(h)
                % set(h(r).label.Parent.Parent,'SizeChangedFcn',[]);
                delete(h.label) 
     
            end
        end
        
    end
    
    
end

function x=roman(n)
% this subfunction converts numbers up to 4999
r = reshape('IVXLCDM   ',2,5);	% the 3 last blank chars are to avoid error for n >= 1000
x = '';
m = floor(log10(n)) + 1;	% m is the number of digit
% n is processed sequentially for each digit
for i = m:-1:1
	ii = fix(n/10^(i-1));	% ii is the digit (0 to 9)
	% Roman numeral is a concatenation of r(1:2,i) and r(1,i+1)
	% combination with regular rules (exception for 4000 = MMMM)
	% Note: the expression uses REPMAT behavior which returns empty
	% string for N <= 0
	x = [x,repmat(r(1,i),1,ii*(ii < 4 | (ii==4 & i==4)) + (ii == 9) + (ii==4 & i < 4)), ...
		   repmat([r(2,i),repmat(r(1,i),1,ii-5)],1,(ii >= 4 & ii <= 8 & i ~= 4)), ...
		   repmat(r(1,i+1),1,(ii == 9))];
	n = n - ii*10^(i-1);	% substract the most significant digit
end
end