function [mdl,int,rts] = fitloglinear(P,c,varargin)

%FITLOGLINEAR fit loglinear model to point pattern
%
% Syntax
%
%     [mdl,int] = fitloglinear(P,c)
%     [mdl,int] = fitloglinear(P,c,pn,pv,...)
%     [mdl,int,mx] = ...
%
% Description
%
%     Loglinear models embrace numerous models that can be fitted to
%     homogeneous and inhomogeneous Poisson processes on river networks. 
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
%                  polynomial of one covariate: 'poly4'
%     
%     In addition, fitloglinear accepts parameter name/value pairs of
%     fitglm or stepwiseglm.
%
% Output arguments
%
%     mdl    model (GeneralizedLinearModel)
%     int    node-attribute list with modelled intensities
%     mx     for higher-order polynomials of a single-variable model, 
%            mx returns the location of maxima in the intensity function. 
%     
% Example: Create an inhomogeneous Poisson process with the intensity being 
%          a function of elevation. Simulate a random point pattern and fit
%          the model.
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
% See also: PPS, PPS/random, fitglm, stepwiseglm
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 15. September, 2021 

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

if nargout == 3
    if mdl.NumPredictors ~= 1
        rts = [];
        return
    end
    
    coeffs = mdl.Coefficients.Estimate(1:end);
    coeffs = flipud(coeffs);
    % first derivative
    p1     = polyder(coeffs);
    % p1     = coeffs.*(numel(coeffs):-1:1)';
    % minima and maxima
    rts    = roots(p1);
    if isempty(rts)
        % there are no minima or maxima
    else
        % second derivative
        p2     = polyder(p1);
        % p2     = p1(1:end-1).*((numel(coeffs)-1):-1:1)';
        y      = polyval(p2,rts);
        rts    = rts(y<0);
    end
end

end

function pnpv = expandstruct(s)

pn = fieldnames(s);
pv = struct2cell(s);

pnpv = [pn pv];
pnpv = pnpv';
pnpv = pnpv(:)';
end