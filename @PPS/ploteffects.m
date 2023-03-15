function ploteffects(P,mdl,varargin)

%PLOTEFFECTS Plot of slices through fitted generalized linear regression
%
% Syntax
%
%     ploteffects(P,mdl)
%     ploteffects(P,mdl,covariate)
%
% Description
%
%     ploteffects plots the individual effects of a loglinear point process
%     model.
%
% Input arguments
%
%     P          PPS object
%     mdl        object of 'GeneralizedLinearModel'
%     covariate  number of covariate to be plotted. For example, if there 
%                two predictor variables in the model mdl, then covariate
%                can be 1 or 2, or [1 2].
%
%     Parameter name value
%
%     'plot'   {true} or false
%     'n'      number of values linear spaced for a covariate {200}
%     'fixedvars'   fixed values (by default, this is the mean of each
%              covariate).
%     'plotintensity' {true} or false. Plots a horizontal line with the
%              intensity of the point pattern
%     'indicators' {false} or true. If true, lines indicating point
%              locations at the bottom of the plot
%     'color'  color of line plots
%     'intcolor'    color of horizontal line with intensity of point
%              pattern (see option plotintensity)
%
% Example: see fitloglinear
%
% 
% See also: PPS, PPS/fitloglinear, PPS/roc 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019
%% Parse inputs
p = inputParser;
p.FunctionName = 'PPS/ploteffects';
% Add elevation
addRequired(p,'P');
addRequired(p,'mdl',@(x) isa(x,'GeneralizedLinearModel'));
addOptional(p,'covariate',1,@(x) numel(x) >= 1 && numel(x) <= (size(mdl.Variables,2)-1));
addParameter(p,'plot',true);
addParameter(p,'n',200);
addParameter(p,'fixedvars',mean(mdl.Variables{:,1:end-1}));
addParameter(p,'plotintensity',true);
addParameter(p,'indicators',false);
addParameter(p,'varnames','');
addParameter(p,'color','k');
addParameter(p,'intcolor',[.4 .4 .4]);
% Parse
parse(p,P,mdl,varargin{:});

covariate = p.Results.covariate;
nvars     = numel(covariate);

% find a nice layout of subplots;
if nvars > 1
nn = numSubplots(nvars);
nrows = nn(1);
ncols = nn(2);
end

fixedvars = p.Results.fixedvars;
fixedvars = repmat(fixedvars,p.Results.n,1);

d   = distance(P.S,'node_to_node');
d   = mean(d);


for r = 1:nvars    
    
    if nvars > 1
    subplot(nrows,ncols,r);
    end
    predictor  = linspace(min(mdl.Variables{:,covariate(r)}),...
                          max(mdl.Variables{:,covariate(r)}),...
                          p.Results.n)';
    preds      = fixedvars;
    preds(:,covariate(r)) = predictor;
    [int,ci] = predict(mdl,preds);
    int = int/d;
    ci  = ci/d;
    
    if p.Results.plot

        plot(predictor,ci,'--','color',p.Results.color);
        hold on
        plot(predictor,int,'-','color',p.Results.color);
        
        if p.Results.plotintensity
            ii = intensity(P);
            plot(xlim,[ii ii],':','color',p.Results.intcolor);
        end
        
        if p.Results.indicators
            xl = mdl.Variables{:,covariate(r)}(P.PP);
            xlinerel(xl,0.03);
        end
        
        hold off
        
        if isempty(p.Results.varnames)  
            lab = ['x_{' num2str(covariate(r)) '}'];
        else
            lab = [p.Results.varnames{r}];
        end
        xlabel(lab);
        ylabel(['\lambda(' lab ')']);
    end
end
    
end


function [p,n]=numSubplots(n)
% function [p,n]=numSubplots(n)
%
% Purpose
% Calculate how many rows and columns of sub-plots are needed to
% neatly display n subplots. 
%
% Inputs
% n - the desired number of subplots.     
%  
% Outputs
% p - a vector length 2 defining the number of rows and number of
%     columns required to show n plots.     
% [ n - the current number of subplots. This output is used only by
%       this function for a recursive call.]
%
%
%
% Example: neatly lay out 13 sub-plots
% >> p=numSubplots(13)
% p = 
%     3   5
% for i=1:13; subplot(p(1),p(2),i), pcolor(rand(10)), end 
%
%
% Rob Campbell - January 2010
   
    
while isprime(n) && n>4
    n=n+1;
end
p=factor(n);
if length(p)==1
    p=[1,p];
    return
end
while length(p)>2
    if length(p)>=4
        p(1)=p(1)*p(end-1);
        p(2)=p(2)*p(end);
        p(end-1:end)=[];
    else
        p(1)=p(1)*p(2);
        p(2)=[];
    end    
    p=sort(p);
end
%Reformat if the column/row ratio is too large: we want a roughly
%square design 
while p(2)/p(1)>2.5
    N=n+1;
    [p,n]=numSubplots(N); %Recursive!
end

end
