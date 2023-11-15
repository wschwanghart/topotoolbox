function clr = ttclrr(name)

%TTCLRR Colors obtained from colornames
%
% Syntax
%
%     clr  = ttclrr(name)
%     dict = ttclrr('all')
%
% Description
%
%     ttclrr returns a three-element rgb color vector with values ranging
%     between 0 and 1. 
%
%     ttclrr('all') returns a dictionary (or table, depending on MATLAB 
%     version) with all colornames and values.
%
%     ttclrr without input and output arguments creates a plot showing all
%     colors.
%
% Input arguments
%
%     name       color name
%
% Output arguments
%
%     clr        1x3 rgb vector (double)
%     dict       dictionary or table with all colors
%
% R-code to generate the csv-file rcolors.csv
%
%     # Extract colors into data.frame
%     
%     # Color names
%     clr <- colors()
%     # Color names to RGB
%     RGB <- col2rgb(clr)
%     # Transpose and divide so that
%     # colors range between 0 and 1
%     RGB <- t(RGB)/255
%     
%     # RGB to data.frame
%     df <- data.frame(clr = clr,
%                      R = RGB[,1],
%                      G = RGB[,2],
%                      B = RGB[,3])
%     write.csv(df,file = "rcolors.csv",
%               row.names = FALSE) 
%
% See also: ttclr, ttscm, ttcmap
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 14. November, 2023 

persistent dict

% Read csv-file with color names and values and make it a persistent
% variable
if isempty(dict)
    p = mfilename("fullpath");
    p = fileparts(p);
    p = fullfile(p,"rcolors.csv");
    clrtable = readtable(p);
    if verLessThan('matlab','9.13')
        dict = clrtable;
    else
        rgb = table2array(clrtable(:,2:end));
        rgb = mat2cell(rgb,ones(size(rgb,1),1),3);
        clrtable.clr = string(clrtable.clr);
        dict = dictionary(clrtable.clr,rgb);
    end
end

if nargin ~= 0

    switch lower(name)
        case 'all'
            clr = dict;
            return
    end

    if verLessThan('matlab','9.13')

        name = validatestring(name,dict.clr);
        I = strcmp(name,dict.clr);
        clr = table2array(dict(I,2:end));

    else

        name = validatestring(name,keys(dict));
        clr = dict(name);
        clr = clr{1};

    end

else

    if verLessThan('matlab','9.13')

        n = size(dict,1);
        names = dict.clr;
        rgb = [dict.R dict.G dict.B];

    else

        n = numEntries(dict);
        names = keys(dict);
        rgb = values(dict);
        rgb = vertcat(rgb{:});

    end

    figure
    h = axes;
    nrows = 50;
    ncols = ceil(n/nrows);
    h.XLim = [1 ncols*2];
    h.YLim = [-1 nrows+1];

    
    
    hold on
    for r = 1:numel(names)
        row = mod(r,nrows);
        col = ceil((r+1)/nrows)*2 - 1;
        y = [row row row-1 row-1 row];
        x = [col col+1 col+1 col col];
        patch(h,x,y,rgb(r,:),'EdgeColor','k');
        text(h,col+1,row,[' ' names{r}],'FontSize',6,'VerticalAlignment','bottom')
    end
    axis ij
    axis off
end

        
