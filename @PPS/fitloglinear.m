function [mdl,int,rts,rtssigma,ismx,sigmapred] = fitloglinear(P,c,varargin)

%FITLOGLINEAR fit loglinear model to point pattern
%
% Syntax
%
%     [mdl,int] = fitloglinear(P,c)
%     [mdl,int] = fitloglinear(P,c,pn,pv,...)
%     [mdl,int,mx,sigmamx,ismx,sigmapred] = ...
%                   fitloglinear(P,c,'modelspec','poly2',...)
%
% Description
%
%     Loglinear models embrace numerous models that can be fitted to
%     homogeneous and inhomogeneous Poisson processes on river networks. 
%     fitloglinear fits a model via logistic regression. Note that this
%     requires that no duplicate points exist. fitloglinear uses fitglm
%     to fit the logistic regression model.
%
% Input arguments
%
%     P      instance of PPS
%     c      covariates (node-attribute lists)
%     
%     Parameter name/value pairs
%     
%     'stepwise'   {false} or true. If true, fitloglinear uses stepwiseglm
%                  to fit the model.
%     'modelspec'  see fitglm. For example, for fitting a forth-order
%                  polynomial: 'poly4'
%     
%     In addition, fitloglinear accepts parameter name/value pairs of
%     fitglm or stepwiseglm.
%
% Output arguments
%
%     mdl          model (GeneralizedLinearModel)
%     int          node-attribute list with modelled intensities
%
%     If (and only if) 'modelspec' is 'poly2', then additional output 
%     arguments are returned
%
%     mx           for second-order polynomials of a single-variable model, 
%                  mx returns the location of maximum in the intensity 
%                  function.
%     sigmamx      the standard error of the location of the maximum
%     ismx         true or false. True, if mx is a maximum, and false
%                  otherwise.
%     sigmapred    the +/- range in which 66% of points will lie, corrected 
%                  for the abundance of the covariate.
%     
% Example 1: Create an inhomogeneous Poisson process with the intensity
%          being a function of elevation. Simulate a random point pattern
%          and fit the model.
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S,1);
%     P = PPS(S,'PP',[],'z',DEM);
%     P = simulate(P,'intensity',getnal(S,DEM)/1e6);
%     subplot(1,2,1)
%     plot(P)
%     subplot(1,2,2)
%     plotdz(P)
%     [mdl,int] = fitloglinear(P,DEM,'model','poly1');
%     figure
%     subplot(1,2,1)
%     plotc(P,int)
%     subplot(1,2,2)
%     ploteffects(P,mdl,1)
%
% Example 2: Model the distribution of knickpoints as a second-order
%            loglinear model.
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S,1);
%     
%     % Find knickpoints
%     [~,kp] = knickpointfinder(S,DEM,'tol',30,...
%        'split',false,...
%        'verbose',false,...
%        'plot',false);
%     P = PPS(S,'PP',kp.IXgrid,'z',DEM);
%     P.PP(~(DEM.Z(S.IXgrid(P.PP)) < 1000)) = [];
%
%     % also try
%     % P.PP(~(DEM.Z(S.IXgrid(P.PP)) < 1000 & ...
%     %        DEM.Z(S.IXgrid(P.PP)) > 800)) = [];
%     
%     A = flowacc(FD);
%     c = chitransform(S,A,'mn',0.4);
%     
%     [mdl,int,mx,sigmamx,ismx,sigmapred] = ...
%                   fitloglinear(P,c,'modelspec','poly2');
%
%     plotdz(P,'distance',c,'ColorData','k')
%     yyaxis right
%     ploteffects(P,mdl,'patch',true)
%     xline([mx-sigmapred mx+sigmapred],':')
%     xline([mx-sigmamx mx+sigmamx],'--')
%     xline(mx)
%     
%
% See also: PPS, PPS/random, fitglm, stepwiseglm
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 9. September, 2022 

p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'stepwise',false,@(x) isscalar(x));
addParameter(p,'modelspec','linear');
addParameter(p,'distribution','binomial');
addParameter(p,'weights',getnal(P.S)+1);
parse(p,varargin{:});

% Get covariates
X = getcovariate(P,c);
% Get response variable
y = zeros(numel(P.S.x),1);
y(P.PP) = 1;

% Options for fitglm and stepwise
glmopts = p.Unmatched;
glmopts = expandstruct(glmopts);

% Handle weights
fllopts = p.Results;
validateattributes(fllopts.weights,{'single','double'},{'>',0,'nonnan','finite'});
if ~isnal(P.S,fllopts.weights)
    % weights should have as many entries as there are points
    if numel(fllopts.weights) ~= npoints(P)
        error('TopoToolbox:fitloglinear','Wrong number of elements in the weights vector')
    end
    w = getnal(P.S)+mean(fllopts.weights);
    w(P.PP) = fllopts.weights;
    fllopts.weights = w;
end


% Covariates can be supplied as table. This needs some handling.
if istable(X)
    tbl = [X table(y,'variablenames',{'response'})];
    inp = {tbl};
else
    inp = {X y};
end


if ~p.Results.stepwise
    
    mdl = fitglm(inp{:},fllopts.modelspec,...
        'Distribution',fllopts.distribution,'weights',fllopts.weights,glmopts{:});
    
else

    % Use stepwise GLM
    mdl = stepwiseglm(inp{:},fllopts.modelspec,...
        'Distribution',fllopts.distribution,'weights',fllopts.weights,glmopts{:});
    
end

% modelled intensities
p   = predict(mdl,X);
d   = distance(P.S,'node_to_node');
d   = mean(d);
int = p./d; %.cellsize;

% --- The following section applies only to second-order loglinear models
%     that have been obtained with 'modelspec','poly2'
if nargout >= 3
    switch lower(fllopts.modelspec) 
        case 'poly2'
        otherwise
            rts = [];
            rtssigma = [];
            ismax = [];
            sigmapred = [];
            warning('TopoToolbox:fitloglinear',...
                ['Third to sixth output only available if\n' ...
                 'model has one variable and modelspec is poly2.'])
        return
    end

    % Coefficient of x
    b1 = mdl.Coefficients.Estimate(2);
	% Coefficient of x^2
    b2 = mdl.Coefficients.Estimate(3);
	% Covariance matrix
    C  = mdl.CoefficientCovariance(2:3,2:3);

    % Location of the maximum
	% b1 + 2*b2*x = 0
    rts = -0.5*b1/b2;
    
	% Standard errors of coefficients
    sb1 = mdl.Coefficients.SE(2);
    sb2 = mdl.Coefficients.SE(3);
	% 
    rtssigma = sqrt(((sb1/b1)^2 + (sb2/b2)^2 - 2*C(2,1)/(b1*b2)))*rts;
    %
    ismx = b2 < 0;

    % Find standard deviation of predictions
    % Can be found be finding the inflection points (zeros of the second
    % derivative of the exponential function
    
    % Let's swap parameter names (so that they are compatible with equations
    % from Wolfram Alpha ;-)
    b3 = b2;
    b2 = b1;
    val1 = (-b2*b3 - sqrt(2)*sqrt(-b3^3))/(2*b3^2);
    % val2 = (sqrt(2)*sqrt(-b3^3) - b2*b3)/(2*b3^2);
    sigmapred = abs(val1-rts);

end
end
    
% ---------------- some helper functions -------------------------------
function pnpv = expandstruct(s)

pn = fieldnames(s);
pv = struct2cell(s);

pnpv = [pn pv];
pnpv = pnpv';
pnpv = pnpv(:)';
end