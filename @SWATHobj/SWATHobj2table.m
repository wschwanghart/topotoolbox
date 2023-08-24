function t = SWATHobj2table(SW,aggfun,varnames)

%SWATHobj2table Convert SWATHobj to table
%
% Syntax
%
%     t = SWATHobj2table(SW)
%     t = SWATHobj2table(SW,'xdist_as_var' or 'ydist_as_var')
%     t = SWATHobj2table(SW,aggfun,varnames)
%
% Description
%
%     SWATHobj2table converts a SWATHobj to a table. 
%     SWATHobj2table(SW) returns a table so that each swath point has one
%     row with the variables X, Y, Z, distx and disty.
%     SWATHobj2table(SW,'xdist_as_var') creates a table in which unique 
%     distances along x ('xdist_as_var') are stored as variable, and
%     z-values are stored as variables according to their y-distance.
%     SWATHobj2table(SW,'ydist_as_var') creates a table in which unique 
%     distances perpendicular to x ('ydist_as_var') are stored as variable, 
%     and z-values are stored as variables according to their distance
%     along the SWATHobj.
%     SWATHobj2table(SW,aggfun,varnames) aggragates z-values perpendicular
%     to the swath distance using the aggregation function(s) in aggfun and
%     stores them in variables with names varnames.
% 
% Input arguments
%
%     SW       SWATHobj
%     aggfun   string of function, anonymous function, or cell array of
%              strings of functions or anonymous functions. Note that it
%              might be necessary to define functions that can handle nans
%              (e.g. nanmean).
%     varnames cell array of variable names (same number as aggfuns)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     A  = flowacc(FD);
%     S = STREAMobj(FD,A>100);
%     S = trunk(klargestconncomps(S,1));
%     [x,y] = STREAMobj2XY(S);
%     ix = ~isnan(x);
%     SW = SWATHobj(DEM,x(ix),y(ix),'smooth',200);
%     
%     t = SWATHobj2table(SW,...
%         {@nanmean @(x)prctile(x,5) @(x)prctile(x,95)},...
%         {'m' 'p05' 'p95'});
%     
%     plot(t.distx,t.m,'k-',t.distx,[t.p05 t.p95],'--')
%     legend('mean','5%-ile', '95%-ile')
%
% See also: SWATHobj, table
%
% Author: Wolfgang Schwanghart 
% Date: 24. August, 2023


if nargin == 1

    % Simple table with all data
    t = table;
    t.X = SW.X(:);
    t.Y = SW.Y(:);
    t.Z = SW.Z(:);
    distx   = repmat(SW.distx,1,numel(SW.disty));
    t.distx = distx(:);
    disty   = repmat(SW.disty,numel(SW.distx),1);
    t.disty = disty(:);
    inan = isnan(t.Z);
    if any(inan)
        t = t(~inan,:);

    end
    return
elseif nargin == 2
    aggfun = validatestring(aggfun,{'ydist_as_var','xdist_as_var'});

    % Table with unique distance values along/perpendicular to swath
    switch aggfun
        case 'ydist_as_var'
            t = table;
            t.disty = SW.disty(:);
            varnames = cellfun(@(x) ['z' num2str(x)],num2cell(SW.distx), ...
                'UniformOutput',false);
            t = [t array2table(SW.Z,"VariableNames",varnames)];
            return
        case 'xdist_as_var'
            t = table;
            t.distx = SW.distx(:);
            varnames = cellfun(@(x) ['z' num2str(x)],num2cell(SW.disty), ...
                'UniformOutput',false);
            t = [t array2table(SW.Z',"VariableNames",varnames)];
            return
    end
end

% aggregating multiple inputs
% place aggregate functions and names in cell array
if ~iscell(aggfun)
    aggfun = {aggfun};
end

if ~iscell(varnames)
    varnames = {varnames};
end

% Do aggfun and varnames have the same number of entries
if length(aggfun) ~= length(varnames)
    error('aggfun and varnames must have the same number of entries.')
end

t = table;
t.distx = SW.distx(:);
for r = 1:numel(aggfun)
    if ischar(aggfun{r}) || isstring(aggfun{r})
        fun = str2func(aggfun{r});

    else
        fun = aggfun{r};
    end

    t.(varnames{r}) = fun(SW.Z)';
end
