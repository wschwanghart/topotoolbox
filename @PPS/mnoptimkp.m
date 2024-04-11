function results = mnoptimkp(P,FD,varargin)

%MNOPTIMKP Find optimal mn ratio based on knickpoint locations
%
% Syntax
%
%     mn = mnoptimkp(P,FD)
%     mn = mnoptimkp(P,FD,'method',method)
%     results = mnoptimkp(P,FD,'method','nlinfit','pn','pv')
%     results = mnoptimkp(P,FD,'method','deviance','pn','pv')
%
% Description
%
%     Knickpoints that originate from a common location cluster in chi
%     space. This function provides access to numerous methods to determine
%     the mn-ratio based on a minimization of the variance of knickpoint
%     chi values. The methods 'nlinfit' and 'deviance' return structure
%     arrays with additional stream-power related parameters and variables.
%
%     The function is currently experimental. This means, that it provides
%     numerous options that are not fully tested nor validated. 
%
% Input arguments
%
%     P       Point pattern object (PPS)
%     FD      FLOWobj
%
%     Parameter name/value pairs (note that some parameters are only valid
%     in specific combinations). Some methods do not take into account that
%     the abundance of knickpoints depends on the topology of the river
%     network. 
%
%     'method'    
%         'var' - minimizes the coefficient of variation of knickpoints 
%                 chi values
%         'robustvar' - minimizes the coefficient of variation of
%                 knickpoint chi values based on a robust estimate of
%                 variance (see function robustvar)
%         'mad0' - minimizes the absolute deviation from the mean
%         'mad1' - minimizes the absolute deviation from the median
%         'nlinfit' - uses a minimization based on nonlinear fitting.
%                 Requires setting 't' and returns both confidence bounds
%                 on mn and K.
%         'deviance' - minimizes the deviance of a second-order loglinear
%                 model (see fitloglinear). If 't' is provided, than the
%                 approach returns an uncertainty estimate of 'K' but not
%                 'mn'. The method will not work if a second-order model
%                 with a maximum cannot be fitted and is most suited if 
%                 there is an obvious and well-defined cluster of
%                 knickpoints in chi space. The approach normalizes for the 
%                 dependence of knickpoint density on network topology.
%
%     'weights' npoints(P)x1 vector of weights for each point in P. Default 
%               is ones(npoints(P),1).
%     'a0'      reference area (default = 1e6 m^2)
%     'mn0'     initial value for mn search. 
%     'optimizemn' {true} or false. If false, then the function will not
%               search for an optimal value of mn but take mn0 (works only
%               for 'method', 'Deviance')
%     't'       time span since knickpoint migration (required for
%               'nlinfit' and 'deviance'. Is not required for finding mn
%               but for finding K and tau. Default is 10000 years.
%     'K0'      initial value for K. Required by 'nlinfit'
%     'inletix' linear index into GRIDobj. Must be located on the stream
%               network P.S and is the location where additional area
%               is added or substracted.
%     'inletA0' initial value for area to be added.
%     'inletfraction' fraction of areas to be added at several locations
%     'Kint'    confidence intervals of K. Default is [0.025 0.975].
%
% Output variables
%
%      results  mn-ratio or structure array. 
%
% Example
%
%      DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%      FD  = FLOWobj(DEM,'preprocess','c');
%      S = STREAMobj(FD,'minarea',1000);
%      S = klargestconncomps(S,1);
%     
%      % Find knickpoints
%      [~,kp] = knickpointfinder(S,DEM,'tol',30,...
%          'split',false,...
%          'verbose',false,...
%          'plot',false);
%      P = PPS(S,'PP',kp.IXgrid,'z',DEM);
%      P.PP(~(DEM.Z(S.IXgrid(P.PP)) < 1000)) = [];
%      
%      results = mnoptimkp(P,FD,'method','deviance',...
%                          't',1e6,'mn0',0.5,...
%                          'sigmat,1e5);
%
%
% See also: PPS, knickpointfinder, STREAMobj/mnoptim, STREAMobj/mnoptimvar
%  
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. April, 2024

p = inputParser;
addParameter(p,'method','robustvar')
addParameter(p,'a0',1);
addParameter(p,'mn0',0.45);
addParameter(p,'optimizemn',true)
addParameter(p,'weights',[],@(x) numel(x) == npoints(P));
addParameter(p,'t',10000);
addParameter(p,'sigmat',0)
addParameter(p,'K0',1e-4);
addParameter(p,'RobustWgtFun',[]);
addParameter(p,'inletix',[]);
addParameter(p,'inletA0',[]);
addParameter(p,'inletfraction',[]);
addParameter(p,'Kint',[0.025 0.975])

parse(p,varargin{:});

A = flowacc(FD);

if isa(A,'GRIDobj')
    a = getnal(P.S,A);
elseif isnal(P.S,A)
    a = A;
else
    error('The second input variable must be a GRIDobj or a node attribute list')
end

% Anonymous function to calculate chi with mn as input
chifun = @(mn) chitransform(P.S,a,'mn',mn,'a0',p.Results.a0);


switch lower(p.Results.method)
    case 'robustvar'
        % Robust variance
        varfun  = @robustcov;
        results = exp(fminsearch(@(mn) chicov(exp(mn)),log(p.Results.mn0))); 
    case 'var'
        % Variance
        varfun  = @var;
        results = exp(fminsearch(@(mn) chicov(exp(mn)),log(p.Results.mn0))); 
    case {'mad0','mad1'}
        v = str2double(p.Results.method(end));
        % mean absolute deviations (mad0)
        % median absolute deviations (mad1)
        results = exp(fminsearch(@(mn) mad(getmarks(P,...
            chifun(exp(mn))/max(chifun(exp(mn)))),v),log(p.Results.mn0)));

    case 'nlinfit'
        % using nlinfit requires t to be known
        % Either use observation weights or RobustWgtFun
        
        % Set options and weights
        if isempty(p.Results.weights)
            options     = statset('RobustWgtFun',p.Results.RobustWgtFun);
            weightinput = cell(0); 
        else
            options     = cell(0);
            weightinput = {'weights',p.Results.weights};
        end
        
        % Function for nlinfit
        nlinfun = @(b,X) getmarks(P,chifun(exp(b(1)))./(exp(b(2)).*(p.Results.a0^exp(b(1)))));
        [beta,R,J,CovB,MSE,ErrorModelInfo] = ...
            nlinfit(zeros(npoints(P),1),zeros(npoints(P),1)+p.Results.t,...
            nlinfun,...
            [log(p.Results.mn0) log(p.Results.K0)],options,...
            weightinput{:});
        
        % Write output to structure
        results.mn = exp(beta(1));
        results.K  = exp(beta(2));
        results.tau = chifun(results.mn) / (results.K * p.Results.a0^results.mn);
        results.taures = R;
        results.J   = J;
        results.CovB = CovB;
        results.MSE  = MSE;
        results.ErrorModelInfo = ErrorModelInfo;
        results.RMSE = sqrt(mean((R-p.Results.t).^2));
        
        
    case 'deviance'
        % using a second-order loglinear model, minimizing deviance of the
        % model
        
        if isempty(p.Results.weights)
            weightinput = cell(0); 
        else
            weightinput = {'weights',p.Results.weights};
        end
        
        if isempty(p.Results.inletix)
            if p.Results.optimizemn
                mn = exp(fminsearch(@(mn) fitloglinear(P,chifun(exp(mn)),...
                    'model','poly2',weightinput{:}).Deviance,...
                    log(p.Results.mn0)));
            else 
                mn = p.Results.mn0;
            end
            c  = chifun(mn);
            [mdl,int,locmax,sigmalocmax] = fitloglinear(P,c,...
                'model','poly2',weightinput{:});
        else
            ix    = p.Results.inletix; % linear index into grid
            [onnetwork,ix]    = ismember(ix,P.S.IXgrid); % ix is now linear index into nal
            if ~all(onnetwork)
                error('Some inlet indices are not located on the stream network.')
            end
            a   = hillslopearea(P.S,FD) + 1; % hillslope area + area of each stream pixel
            add = getnal(P.S);
            if ~isempty(p.Results.inletfraction)
                add(ix) = p.Results.inletfraction;
            else
                add(ix) = ones(numel(ix),1)/numel(ix);
            end
            chifun = @(mn,addA) chitransform(P.S,cumsum(P.S,a + add*addA,'downstream'),'mn',mn,'a0',p.Results.a0);
                
            if p.Results.optimizemn
                mna = fminsearchbnd(@(mna) fitloglinear(P,chifun(exp(mna(1)),mna(2)),...
                    'model','poly2',weightinput{:}).Deviance,...
                    [log(p.Results.mn0) (p.Results.inletA0/(P.S.cellsize^2))]);
                mna(1) = exp(mna(1));
            else
                mn  = p.Results.mn0;
                mna = fminsearch(@(mna) fitloglinear(P,chifun(mn,mna),...
                    'model','poly2',weightinput{:}).Deviance,...
                    (p.Results.inletA0/(P.S.cellsize^2)));
                mna = [mn mna];
            end
            c      = chifun(mna(1),mna(2));
            [mdl,int,locmax,sigmalocmax] = fitloglinear(P,c,...
                'model','poly2',weightinput{:});
            mn = mna(1);
            addA = mna(2);
            
        end
        
        results.mn  = mn;
        results.chi = c;
        results.chifun = chifun;
        if ~isempty(p.Results.inletix)
            results.addA_pix = addA;
            results.addA_m2  = addA*P.S.cellsize^2;
            results.addA_km2 = results.addA_m2 / 1e6;
        end
        results.mdl = mdl;
        results.int = int;
        results.chimax = locmax;
        
        if ~isempty(p.Results.t)
            tini = p.Results.t; % onset of incision
            stini = p.Results.sigmat; % standard deviation of the onset
            a0    = p.Results.a0;
            mn    = results.mn;

            results.K      = results.chimax/(tini * a0^mn);
            if stini == 0
                % Without error of incision
                results.Ksigma = sigmalocmax/(tini * a0^mn); 
            else
                % With error of incision
                results.Ksigma = hypot(stini/tini,sigmalocmax/results.chimax) * (sigmalocmax/(tini*a0^mn));
            end
            cdfvals      = p.Results.Kint;
            results.Kint = results.Ksigma*icdf('normal',cdfvals,0,1) + results.K;
            % results.Kint = (results.chimax + results.Ksigma*icdf('normal',cdfvals,0,1))/(tini *a0^mn); 
            results.tau = results.chi / (results.K * a0^mn);
            results.tauint = results.chi ./ (results.Kint *a0^mn);
        end
        
end

function cofv = chicov(mn)
% Coefficient of variation
c = chifun(mn);
% normalize chi to range between 0 and 1
cn = c/max(c);
cofv = varfun(getmarks(P,cn));
%cofv = varfun(getmarks(P,c))/mean(getmarks(P,c));

end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to get error ellipse from covariance matrix
function h=error_ellipse(varargin)
% ERROR_ELLIPSE - plot an error ellipse, or ellipsoid, defining confidence region
%    ERROR_ELLIPSE(C22) - Given a 2x2 covariance matrix, plot the
%    associated error ellipse, at the origin. It returns a graphics handle
%    of the ellipse that was drawn.
%
%    ERROR_ELLIPSE(C33) - Given a 3x3 covariance matrix, plot the
%    associated error ellipsoid, at the origin, as well as its projections
%    onto the three axes. Returns a vector of 4 graphics handles, for the
%    three ellipses (in the X-Y, Y-Z, and Z-X planes, respectively) and for
%    the ellipsoid.
%
%    ERROR_ELLIPSE(C,MU) - Plot the ellipse, or ellipsoid, centered at MU,
%    a vector whose length should match that of C (which is 2x2 or 3x3).
%
%    ERROR_ELLIPSE(...,'Property1',Value1,'Name2',Value2,...) sets the
%    values of specified properties, including:
%      'C' - Alternate method of specifying the covariance matrix
%      'mu' - Alternate method of specifying the ellipse (-oid) center
%      'conf' - A value betwen 0 and 1 specifying the confidence interval.
%        the default is 0.5 which is the 50% error ellipse.
%      'scale' - Allow the plot the be scaled to difference units.
%      'style' - A plotting style used to format ellipses.
%      'clip' - specifies a clipping radius. Portions of the ellipse, -oid,
%        outside the radius will not be shown.
%
%    NOTES: C must be positive definite for this function to work properly.
default_properties = struct(...
  'C', [], ... % The covaraince matrix (required)
  'mu', [], ... % Center of ellipse (optional)
  'conf', 0.5, ... % Percent confidence/100
  'scale', 1, ... % Scale factor, e.g. 1e-3 to plot m as km
  'style', '', ...  % Plot style
  'clip', inf); % Clipping radius
if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.C = varargin{1};
  varargin(1) = [];
end
if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.mu = varargin{1};
  varargin(1) = [];
end
if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.conf = varargin{1};
  varargin(1) = [];
end
if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.scale = varargin{1};
  varargin(1) = [];
end
if length(varargin) >= 1 & ~ischar(varargin{1})
  error('Invalid parameter/value pair arguments.') 
end
prop = getopt(default_properties, varargin{:});
C = prop.C;
if isempty(prop.mu)
  mu = zeros(length(C),1);
else
  mu = prop.mu;
end
conf = prop.conf;
scale = prop.scale;
style = prop.style;
if conf <= 0 | conf >= 1
  error('conf parameter must be in range 0 to 1, exclusive')
end
[r,c] = size(C);
if r ~= c | (r ~= 2 & r ~= 3)
  error(['Don''t know what to do with ',num2str(r),'x',num2str(c),' matrix'])
end
x0=mu(1);
y0=mu(2);
% Compute quantile for the desired percentile
k = sqrt(qchisq(conf,r)); % r is the number of dimensions (degrees of freedom)
hold_state = get(gca,'nextplot');
if r==3 & c==3
  z0=mu(3);
  
  % Make the matrix has positive eigenvalues - else it's not a valid covariance matrix!
  if any(eig(C) <=0)
    error('The covariance matrix must be positive definite (it has non-positive eigenvalues)')
  end
  % C is 3x3; extract the 2x2 matricies, and plot the associated error
  % ellipses. They are drawn in space, around the ellipsoid; it may be
  % preferable to draw them on the axes.
  Cxy = C(1:2,1:2);
  Cyz = C(2:3,2:3);
  Czx = C([3 1],[3 1]);
  [x,y,z] = getpoints(Cxy,prop.clip);
  h1=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
  [y,z,x] = getpoints(Cyz,prop.clip);
  h2=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
  [z,x,y] = getpoints(Czx,prop.clip);
  h3=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
  
  [eigvec,eigval] = eig(C);
  [X,Y,Z] = ellipsoid(0,0,0,1,1,1);
  XYZ = [X(:),Y(:),Z(:)]*sqrt(eigval)*eigvec';
  
  X(:) = scale*(k*XYZ(:,1)+x0);
  Y(:) = scale*(k*XYZ(:,2)+y0);
  Z(:) = scale*(k*XYZ(:,3)+z0);
  h4=surf(X,Y,Z);
  colormap gray
  alpha(0.3)
  camlight
  if nargout
    h=[h1 h2 h3 h4];
  end
elseif r==2 & c==2
  % Make the matrix has positive eigenvalues - else it's not a valid covariance matrix!
  if any(eig(C) <=0)
    error('The covariance matrix must be positive definite (it has non-positive eigenvalues)')
  end
  [x,y,z] = getpoints(C,prop.clip);
  h1=plot(scale*(x0+k*x),scale*(y0+k*y),prop.style);
  set(h1,'zdata',z+1)
  if nargout
    h=h1;
  end
else
  error('C (covaraince matrix) must be specified as a 2x2 or 3x3 matrix)')
end
%axis equal
set(gca,'nextplot',hold_state);
end
%---------------------------------------------------------------
% getpoints - Generate x and y points that define an ellipse, given a 2x2
%   covariance matrix, C. z, if requested, is all zeros with same shape as
%   x and y.
function [x,y,z] = getpoints(C,clipping_radius)
n=100; % Number of points around ellipse
p=0:pi/n:2*pi; % angles around a circle
[eigvec,eigval] = eig(C); % Compute eigen-stuff
xy = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
x = xy(:,1);
y = xy(:,2);
z = zeros(size(x));
% Clip data to a bounding radius
if nargin >= 2
  r = sqrt(sum(xy.^2,2)); % Euclidian distance (distance from center)
  x(r > clipping_radius) = nan;
  y(r > clipping_radius) = nan;
  z(r > clipping_radius) = nan;
end
end
%---------------------------------------------------------------
function x=qchisq(P,n)
% QCHISQ(P,N) - quantile of the chi-square distribution.
if nargin<2
  n=1;
end
s0 = P==0;
s1 = P==1;
s = P>0 & P<1;
x = 0.5*ones(size(P));
x(s0) = -inf;
x(s1) = inf;
x(~(s0|s1|s))=nan;
for ii=1:14
  dx = -(pchisq(x(s),n)-P(s))./dchisq(x(s),n);
  x(s) = x(s)+dx;
  if all(abs(dx) < 1e-6)
    break;
  end
end
end
%---------------------------------------------------------------
function F=pchisq(x,n)
% PCHISQ(X,N) - Probability function of the chi-square distribution.
if nargin<2
  n=1;
end
F=zeros(size(x));
if rem(n,2) == 0
  s = x>0;
  k = 0;
  for jj = 0:n/2-1;
    k = k + (x(s)/2).^jj/factorial(jj);
  end
  F(s) = 1-exp(-x(s)/2).*k;
else
  for ii=1:numel(x)
    if x(ii) > 0
      F(ii) = quadl(@dchisq,0,x(ii),1e-6,0,n);
    else
      F(ii) = 0;
    end
  end
end
end
%---------------------------------------------------------------
function f=dchisq(x,n)
% DCHISQ(X,N) - Density function of the chi-square distribution.
if nargin<2
  n=1;
end
f=zeros(size(x));
s = x>=0;
f(s) = x(s).^(n/2-1).*exp(-x(s)/2)./(2^(n/2)*gamma(n/2));
end
%---------------------------------------------------------------
function properties = getopt(properties,varargin)
%GETOPT - Process paired optional arguments as 'prop1',val1,'prop2',val2,...
%
%   getopt(properties,varargin) returns a modified properties structure,
%   given an initial properties structure, and a list of paired arguments.
%   Each argumnet pair should be of the form property_name,val where
%   property_name is the name of one of the field in properties, and val is
%   the value to be assigned to that structure field.
%
%   No validation of the values is performed.
%
% EXAMPLE:
%   properties = struct('zoom',1.0,'aspect',1.0,'gamma',1.0,'file',[],'bg',[]);
%   properties = getopt(properties,'aspect',0.76,'file','mydata.dat')
% would return:
%   properties = 
%         zoom: 1
%       aspect: 0.7600
%        gamma: 1
%         file: 'mydata.dat'
%           bg: []
%
% Typical usage in a function:
%   properties = getopt(properties,varargin{:})
% Process the properties (optional input arguments)
prop_names = fieldnames(properties);
TargetField = [];
for ii=1:length(varargin)
  arg = varargin{ii};
  if isempty(TargetField)
    if ~ischar(arg)
      error('Propery names must be character strings');
    end
    f = find(strcmp(prop_names, arg));
    if length(f) == 0
      error('%s ',['invalid property ''',arg,'''; must be one of:'],prop_names{:});
    end
    TargetField = arg;
  else
    % properties.(TargetField) = arg; % Ver 6.5 and later only
    properties = setfield(properties, TargetField, arg); % Ver 6.1 friendly
    TargetField = '';
  end
end
if ~isempty(TargetField)
  error('Property names and values must be specified in pairs.');
end

end

%% 
function [x,fval,exitflag,output] = fminsearchbnd(fun,x0,LB,UB,options,varargin)
% FMINSEARCHBND: FMINSEARCH, but with bound constraints by transformation
% usage: x=FMINSEARCHBND(fun,x0)
% usage: x=FMINSEARCHBND(fun,x0,LB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options,p1,p2,...)
% usage: [x,fval,exitflag,output]=FMINSEARCHBND(fun,x0,...)
% 
% arguments:
%  fun, x0, options - see the help for FMINSEARCH
%
%  LB - lower bound vector or array, must be the same size as x0
%
%       If no lower bounds exist for one of the variables, then
%       supply -inf for that variable.
%
%       If no lower bounds at all, then LB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
%  UB - upper bound vector or array, must be the same size as x0
%
%       If no upper bounds exist for one of the variables, then
%       supply +inf for that variable.
%
%       If no upper bounds at all, then UB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
% Notes:
%
%  If options is supplied, then TolX will apply to the transformed
%  variables. All other FMINSEARCH parameters should be unaffected.
%
%  Variables which are constrained by both a lower and an upper
%  bound will use a sin transformation. Those constrained by
%  only a lower or an upper bound will use a quadratic
%  transformation, and unconstrained variables will be left alone.
%
%  Variables may be fixed by setting their respective bounds equal.
%  In this case, the problem will be reduced in size for FMINSEARCH.
%
%  The bounds are inclusive inequalities, which admit the
%  boundary values themselves, but will not permit ANY function
%  evaluations outside the bounds. These constraints are strictly
%  followed.
%
%  If your problem has an EXCLUSIVE (strict) constraint which will
%  not admit evaluation at the bound itself, then you must provide
%  a slightly offset bound. An example of this is a function which
%  contains the log of one of its parameters. If you constrain the
%  variable to have a lower bound of zero, then FMINSEARCHBND may
%  try to evaluate the function exactly at zero.
%
%
% Example usage:
% rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
%
% fminsearch(rosen,[3 3])     % unconstrained
% ans =
%    1.0000    1.0000
%
% fminsearchbnd(rosen,[3 3],[2 2],[])     % constrained
% ans =
%    2.0000    4.0000
%
% See test_main.m for other examples of use.
%
%
% See also: fminsearch, fminspleas
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 4
% Release date: 7/23/06
% size checks
xsize = size(x0);
x0 = x0(:);
n=length(x0);
if (nargin<3) || isempty(LB)
  LB = repmat(-inf,n,1);
else
  LB = LB(:);
end
if (nargin<4) || isempty(UB)
  UB = repmat(inf,n,1);
else
  UB = UB(:);
end
if (n~=length(LB)) || (n~=length(UB))
  error 'x0 is incompatible in size with either LB or UB.'
end
% set default options if necessary
if (nargin<5) || isempty(options)
  options = optimset('fminsearch');
end
% stuff into a struct to pass around
params.args = varargin;
params.LB = LB;
params.UB = UB;
params.fun = fun;
params.n = n;
% note that the number of parameters may actually vary if 
% a user has chosen to fix one or more parameters
params.xsize = xsize;
params.OutputFcn = [];
% 0 --> unconstrained variable
% 1 --> lower bound only
% 2 --> upper bound only
% 3 --> dual finite bounds
% 4 --> fixed variable
params.BoundClass = zeros(n,1);
for i=1:n
  k = isfinite(LB(i)) + 2*isfinite(UB(i));
  params.BoundClass(i) = k;
  if (k==3) && (LB(i)==UB(i))
    params.BoundClass(i) = 4;
  end
end
% transform starting values into their unconstrained
% surrogates. Check for infeasible starting guesses.
x0u = x0;
k=1;
for i = 1:n
  switch params.BoundClass(i)
    case 1
      % lower bound only
      if x0(i)<=LB(i)
        % infeasible starting value. Use bound.
        x0u(k) = 0;
      else
        x0u(k) = sqrt(x0(i) - LB(i));
      end
      
      % increment k
      k=k+1;
    case 2
      % upper bound only
      if x0(i)>=UB(i)
        % infeasible starting value. use bound.
        x0u(k) = 0;
      else
        x0u(k) = sqrt(UB(i) - x0(i));
      end
      
      % increment k
      k=k+1;
    case 3
      % lower and upper bounds
      if x0(i)<=LB(i)
        % infeasible starting value
        x0u(k) = -pi/2;
      elseif x0(i)>=UB(i)
        % infeasible starting value
        x0u(k) = pi/2;
      else
        x0u(k) = 2*(x0(i) - LB(i))/(UB(i)-LB(i)) - 1;
        % shift by 2*pi to avoid problems at zero in fminsearch
        % otherwise, the initial simplex is vanishingly small
        x0u(k) = 2*pi+asin(max(-1,min(1,x0u(k))));
      end
      
      % increment k
      k=k+1;
    case 0
      % unconstrained variable. x0u(i) is set.
      x0u(k) = x0(i);
      
      % increment k
      k=k+1;
    case 4
      % fixed variable. drop it before fminsearch sees it.
      % k is not incremented for this variable.
  end
  
end
% if any of the unknowns were fixed, then we need to shorten
% x0u now.
if k<=n
  x0u(k:n) = [];
end
% were all the variables fixed?
if isempty(x0u)
  % All variables were fixed. quit immediately, setting the
  % appropriate parameters, then return.
  
  % undo the variable transformations into the original space
  x = xtransform(x0u,params);
  
  % final reshape
  x = reshape(x,xsize);
  
  % stuff fval with the final value
  fval = feval(params.fun,x,params.args{:});
  
  % fminsearchbnd was not called
  exitflag = 0;
  
  output.iterations = 0;
  output.funcCount = 1;
  output.algorithm = 'fminsearch';
  output.message = 'All variables were held fixed by the applied bounds';
  
  % return with no call at all to fminsearch
  return
end
% Check for an outputfcn. If there is any, then substitute my
% own wrapper function.
if ~isempty(options.OutputFcn)
  params.OutputFcn = options.OutputFcn;
  options.OutputFcn = @outfun_wrapper;
end
% now we can call fminsearch, but with our own
% intra-objective function.
[xu,fval,exitflag,output] = fminsearch(@intrafun,x0u,options,params);
% undo the variable transformations into the original space
x = xtransform(xu,params);
% final reshape to make sure the result has the proper shape
x = reshape(x,xsize);
% Use a nested function as the OutputFcn wrapper
  function stop = outfun_wrapper(x,varargin);
    % we need to transform x first
    xtrans = xtransform(x,params);
    
    % then call the user supplied OutputFcn
    stop = params.OutputFcn(xtrans,varargin{1:(end-1)});
    
  end
end % mainline end
% ======================================
% ========= begin subfunctions =========
% ======================================
function fval = intrafun(x,params)
% transform variables, then call original function
% transform
xtrans = xtransform(x,params);
% and call fun
fval = feval(params.fun,reshape(xtrans,params.xsize),params.args{:});
end % sub function intrafun end
% ======================================
function xtrans = xtransform(x,params)
% converts unconstrained variables into their original domains
xtrans = zeros(params.xsize);
% k allows some variables to be fixed, thus dropped from the
% optimization.
k=1;
for i = 1:params.n
  switch params.BoundClass(i)
    case 1
      % lower bound only
      xtrans(i) = params.LB(i) + x(k).^2;
      
      k=k+1;
    case 2
      % upper bound only
      xtrans(i) = params.UB(i) - x(k).^2;
      
      k=k+1;
    case 3
      % lower and upper bounds
      xtrans(i) = (sin(x(k))+1)/2;
      xtrans(i) = xtrans(i)*(params.UB(i) - params.LB(i)) + params.LB(i);
      % just in case of any floating point problems
      xtrans(i) = max(params.LB(i),min(params.UB(i),xtrans(i)));
      
      k=k+1;
    case 4
      % fixed variable, bounds are equal, set it at either bound
      xtrans(i) = params.LB(i);
    case 0
      % unconstrained variable.
      xtrans(i) = x(k);
      
      k=k+1;
  end
end
end % sub function xtransform end

