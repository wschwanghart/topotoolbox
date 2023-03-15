function s = rhohat(P,varargin)

%RHOHAT nonparametric estimation of point pattern dependence on covariate
%
% Syntax
%
%     rhohat(P)
%     rhohat(P,pn,pv,...)
%     r = rhohat(...)
%
% Description
%
%     rhohat estimates the dependence of a spatial process on a spatial
%     covariate using a nonparametric approach. The approach is similar to
%     kernel density estimation yet differs in that it accounts for the
%     network structure (the stream network) of the observation window of 
%     the point pattern. 
%
% Input arguments
%
%     P     instance of PPS
%
%     Parameter name/value pairs
%
%     'covariate'    GRIDobj or node-attribute list of a covariate. By
%                    default, the covariate is the distance from the outlet
%                    (P.S.distance)
%     'intline'      {false} or true. If true, then a horizontal line with 
%                    rho = intensity(P) will be plotted.
%     'nsim'         {1000}. Number of bootstrap samples
%     'confintervals'  {true} or false. If true, bootstrapped onfidence 
%                    intervals will be calculated and plotted.
%     'alpha'        {0.05}. Width of the confidence intervals (1-alpha)
%     'ksdensity'    {true} or false. Use a gaussian kernel density to 
%                    estimate rho (uses the function ksdensity)
%     'bandwidth'    bandwidth of the gaussian kernel in ksdensity.
%     'weights'      vector with weights for each point. Values must be
%                    positive.
%     'name'         variable name. Default is 'x'.
%     'indicators'   {true} or false. If true, point indicators will be
%                    plotted on the x-axis.
%     'indlocation'  {'bottom'} or 'top'. Indicator position at the bottom
%                    or at the top. 
%     'indcolor'     {'b'}. Indicator color.
%     'FaceColor'    Color of the confidence bounds.
%     'FaceAlpha'    Transparency of the confidence bounds.
%
% Output arguments
%
%     r             structure array with following fields
%     .rho          estimated intensity
%     .rhol, rhou   lower and upper bootstrapped confidence intervals
%                   around estimated intensities
%     .covariate    covariate values 
%     .rhopp, .rholpp, .rhoupp 
%                   estimated intensities and upper and lower bounds at
%                   points
%     .bandwidth    bandwidth
%  
% Example: Calculate nonparametric dependence of knickpoint locations on
%          chi
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S,1);
%     [~,kp] = knickpointfinder(S,DEM,'tol',30,'split',false);
%     A = flowacc(FD);
%     c = chitransform(S,A);
%     P = PPS(S,'PP',kp.IXgrid,'z',DEM);
%     r = rhohat(P,'cov',c,'name','\chi');
%
% Reference 
%
%     Baddeley A, Chang Y-M, Song Y, Turner R. 2012. Nonparametric
%     estimation of the dependence of a spatial point process on
%     spatial covariates. Statistics and Its Interface 5 : 221–236.
%     DOI: 10.4310/SII.2012.v5.n2.a7
% 
% See also: PPS 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 14. September, 2021

% Check input arguments
p = inputParser;
p.FunctionName = 'rhohat';
addRequired(p,'P',@(x) isa(x,'PPS'));

% Standard color
clr = [0.5059 0.8471 0.8157];

% Add elevation
addParameter(p,'confintervals',true);
addParameter(p,'intline',false);
addParameter(p,'nsim',1000);
addParameter(p,'covariate',P.S.distance);
addParameter(p,'plot',true);
addParameter(p,'alpha',0.05);
addParameter(p,'ksdensity',true);
addParameter(p,'FaceColor',clr);
addParameter(p,'FaceAlpha',0.5);
addParameter(p,'Color','k')
addParameter(p,'LineStyle','-');
addParameter(p,'LineWidth',1.5);
addParameter(p,'weights',ones(npoints(P),1),@(x) numel(x) == npoints(P) & all(x > 0));
addParameter(p,'name','x');
addParameter(p,'indicators',true);
addParameter(p,'indlocation','bottom',@(x) ischar(validatestring(x,{'bottom','top'})));
addParameter(p,'indcolor','b');
addParameter(p,'indlength',0.03);
addParameter(p,'bandwidth',[]);

% Parse
parse(p,P,varargin{:});

c = getcovariate(P,p.Results.covariate);

[N,x,Nb,Nu,Nl,bw] = ecdf(P,'nsim',p.Results.nsim,'covariate',c,...
    'alpha',p.Results.alpha,'bandwidth',p.Results.bandwidth,...
    'confintervals',p.Results.confintervals,...
    'accintervals',false,'plot',false,'ksdensity',p.Results.ksdensity,...
    'function','pdf',...
    'weights',p.Results.weights);

in = intensity(P);
rho  = in*N./Nb;
rhou = in*Nu./Nb;
rhol = in*Nl./Nb;

if p.Results.plot
    
    ax = gca;
    hold(ax,'on');
    box on
    
    % plot intensity line with acceptance intervals
    if p.Results.intline
        [~,xacc,Nbacc,Nuacc,Nlacc] = ecdf(P,'nsim',p.Results.nsim,'covariate',c,...
            'alpha',p.Results.alpha,'bandwidth',p.Results.bandwidth,...
            'confintervals',false,...
            'accintervals',true,'plot',false,'ksdensity',p.Results.ksdensity,'function','pdf');
        rhoaccu = in*Nuacc./Nbacc;
        rhoaccl = in*Nlacc./Nbacc;
        patch([xacc;flipud(xacc)],[rhoaccu;flipud(rhoaccl)],[0.9 0.9 0.9],...
            'EdgeColor','none','FaceAlpha',p.Results.FaceAlpha);
        plot([min(x),max(x)],[in in],'--k');
    end
    
    % plot confidence intervals
    if p.Results.confintervals
        patch([x; flipud(x)],[rhou;flipud(rhol)],p.Results.FaceColor,...
            'EdgeColor','none','FaceAlpha',p.Results.FaceAlpha);
    end
    
    % plot rho line
    plot(x,rho(:),'LineWidth',p.Results.LineWidth,...
                  'LineStyle',p.Results.LineStyle,...
                  'Color',p.Results.Color,...
                  'Marker','none')
    
    % plot indicators
    if p.Results.indicators
        indloc = validatestring(p.Results.indlocation,{'bottom','top'});
        switch indloc
            case 'bottom'
                rel = p.Results.indlength;
            case 'top'
                rel = [1-p.Results.indlength 1];
        end
        xlinerel(c(P.PP),rel,'-','Color',p.Results.indcolor)
        
    end
    hold off
    ylabel(['\rho(' p.Results.name ')'])
    xlabel(p.Results.name)
    
end

if nargout > 0
    s.covariate = c;
    s.rho  = interp1(x,rho,c);
    s.rhol = interp1(x,rhol,c);
    s.rhou = interp1(x,rhou,c);
    s.rhopp = s.rho(P.PP);
    s.rholpp = s.rhol(P.PP);
    s.rhoupp = s.rhou(P.PP);
    s.bandwidth = bw;
    
end