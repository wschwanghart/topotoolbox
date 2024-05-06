function [zs,output] = quantcarve(S,DEM,tau,varargin)

%QUANTCARVE quantile carving
%
% Syntax
%
%     zs = quantcarve(S,DEM,tau)
%     zs = quantcarve(S,z,tau)
%     zs = quantcarve(...,pn,pv,...)
%     [zs,output] = ...
%
%
% Description
%
%     Elevation values along stream networks are frequently affected by
%     large scatter, often as a result of data artifacts or errors. This
%     function returns a node attribute list of elevations calculated by
%     carving the DEM. Conversely to conventional carving, quantcarve will
%     not run along minimas of the DEM. Instead, quantcarve returns a
%     profile that runs along the tau's quantile of elevation conditional 
%     horizontal distance of the river profile.
%
%     The function uses linprog from the Optimization Toolbox.
%
% Input parameters
%
%     S      STREAMobj
%     DEM    Digital elevation model (GRIDobj)
%     tau    quantile (default is 0.5)
%     
%     Parameter name/value pairs {default}
%
%     'mingradient'   positive scalar {0}. Minimum downward gradient. 
%                     Choose carefully, because length profile may dip to
%                     steeply. 
%     'split'         {true} or false. If set to true, quantcarve will
%                     split the network into individual drainage basins and 
%                     process them in parallel.
%     'waitbar'       {true} or false. Applies only if split option is set
%                     to 2. 
%
%
% Output parameters
%
%     zs       node attribute list with smoothed elevation values
%     output   structure array with information about the optimization
%              progress
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     S = trunk(S);
%     zs50  = quantcarve(S,DEM,0.5);
%     zs90  = quantcarve(S,DEM,0.9);
%     zs10  = quantcarve(S,DEM,0.1);
%     plotdz(S,DEM,'color',[0.6 0.6 0.6])
%     hold on
%     plotdzshaded(S,[zs90 zs10]);
%     plotdz(S,zs50,'color','k','LineWidth',1.5)
%     hold off
%
% Algorithm
%
%     This algorithm uses quantile carving to smooth the data. The
%     algorithm is described in Schwanghart and Scherler (2017) (Eq. A12).
%     
% References
%
%     Schwanghart, W., Scherler, D., 2017. Bumps in river profiles: 
%     uncertainty assessment and smoothing using quantile regression 
%     techniques. Earth Surface Dynamics, 5, 821-839. 
%     [DOI: 10.5194/esurf-5-821-2017]
%
% See also: STREAMobj/mincosthydrocon, STREAMobj/crs, STREAMobj/imposemin
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 25. September, 2017

% check and parse inputs
narginchk(2,inf)

if nargin == 2
    tau = 0.5;
end

p = inputParser;
p.FunctionName = 'STREAMobj/quantcarve';
addParameter(p,'split',2);
addParameter(p,'mingradient',0,@(x) isscalar(x) && x>=0);
addParameter(p,'fixedoutlet',false);
addParameter(p,'waitbar',true);
parse(p,varargin{:});

validateattributes(tau,{'numeric'},{'>',0,'<',1},'STREAMobj/quantcarve','tau',3);

% get node attribute list with elevation values
if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    z = getnal(S,DEM);
elseif isnal(S,DEM)
    z = DEM;
else
    error('Imcompatible format of second input argument')
end

if any(isnan(z))
    error('DEM or z may not contain any NaNs')
end

z = double(z);

%% Run in parallel
if p.Results.split == 1
    % Option 1: Process each drainage basin independently
    params = p.Results;
    params.split = false;
    [CS,locS] = STREAMobj2cell(S);
    if numel(CS) > 1
        % run only in parallel if more than one drainage basin
        Cz = cellfun(@(ix) z(ix),locS,'UniformOutput',false);
        Czs = cell(size(CS));
        parfor r = 1:numel(CS)
            Czs{r} = quantcarve(CS{r},Cz{r},tau,params);
        end
        
        zs = nan(size(z));
        for r = 1:numel(CS)
            zs(locS{r}) = Czs{r};
        end
        return
    else
        zs = quantcarve(S,z,tau,params);
        return
    end
    
elseif p.Results.split == 2
    % Option 2: Process tributaries in parallel
    
    [CS,locb,CID] = STREAMobj2cell(S,'trib');
    
    params = p.Results;
    params.split = 0;
    params.fixedoutlet = false;
    wb     = params.waitbar;
    params.waitbar = false;
    
    if wb
        h = waitbar(0,'Processing');
    end
    
    
    ntribs = max(CID);
    
    for r = 1:ntribs
        if wb
            waitbar(r/ntribs,h);
        end
        ii = CID == r;
        if r > 1
            params.fixedoutlet = true;
        end
        CStemp   = CS(ii);
        locbtemp = locb(ii);
        ztribs    = cellfun(@(ix) z(ix),locbtemp,'UniformOutput',false);
        Czstemp   = cell(numel(ii),1);
        
        parfor r2 = 1:numel(CStemp)
            % quantcarve(Stribs{r},ztribs{r},tau,params);
            Czstemp{r2} = quantcarve(CStemp{r2},ztribs{r2},tau,params);
        end
        for r2 = 1:numel(CStemp)
            z(locbtemp{r2}) = Czstemp{r2};
        end
    end
    
    if wb
        close(h)
    end
    
    zs = z;
    return
    
end




%% Carve function starts here
% upstream distance
d  = S.distance;
% nr of nodes
n  = numel(S.IXgrid);

f   = [tau*ones(n,1);(1-tau)*ones(n,1);zeros(n,1)];
% Equalities
if ~p.Results.fixedoutlet
    Aeq = [speye(n),-speye(n),speye(n)];
else 
    OUTL = streampoi(S,'outlet','logical');
    P    = spdiags(+(~OUTL),0,n,n);
    Aeq  = [P,-P,speye(n)];
end

beq = z;
lb  = [zeros(n,1);zeros(n,1);-inf*ones(n,1)];


% gradient constraint
d = 1./(d(S.ix)-d(S.ixc));
A = [sparse(n,n*2) (sparse(S.ix,S.ixc,d,n,n)-sparse(S.ix,S.ix,d,n,n))];

if p.Results.mingradient~=0
    b = zeros(n,1);
    b(S.ix) = -p.Results.mingradient;
else
    b = sparse(n,1);
end


%% Solve the linear programme
% set options
options = optimset('Display','off','algorithm','interior-point'); %'OptimalityTolerance',1e-6,

[bhat,~,~,output] = linprog(f,A,b,Aeq,beq,lb,[],options);
zs = bhat(2*n+1:end);


