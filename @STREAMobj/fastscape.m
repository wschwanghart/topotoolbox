function z = fastscape(S,z,a,varargin)

%FASTSCAPE  Simulation of the stream power incision model
%
% Syntax
%
%     z = fastscape(S,DEM,A)
%     z = fastscape(S,DEM,A,'pn',pv,...)
%
% Description
%
%     1D-simulation of the stream power incision model using the implicit
%     solver commonly known as fastscape. The solver was developed by 
%     Hergarten (2002), Hergarten and Neugebauer (2001) and Braun and 
%     Willett (2013).
%
% Input arguments
%
%     S           STREAMobj 
%     DEM         node-attribute list (nal) or GRIDobj with initial 
%                 elevation values [m]                 
%     A           nal or GRIDobj with upstream areas as returned by flowacc
%                 (fastscape converts pixel counts to square meters, see
%                 also option convertpx
%
%     Parameter name/value pairs
%
%     u           scalar, nal or GRIDobj with uplift rates [m/y]
%     k           sclaar, nal or GRIDobj with erosivities [m^(1-2m) y^(-1)]
%     convertpx   fastscape by default assumes that drainage areas in
%                 a are provided as pixels. Set correctcellsize to false if 
%                 you don't want that the values in a are converted to 
%                 metric units (default = true).
%     bc          baselevel change rate [m/y]. Set negative values if you
%                 want baselevel to drop. Baselevel increase is not
%                 foreseen unless they are equal or less than uplift rate.
%                 By default baselevel is held constant at the elevation of
%                 the outlet(s). You can also define time-variable boundary
%                 conditions using a table with a column t and dz or z. See
%                 the examples below for the syntax.
%     bctype      'rate' or 'elev'
%     m           area exponent (default is 0.5)
%     n           slope exponent (default is 1)
%     tspan       simulation time [y] (default is 100000)
%     dt          time step length [y] (default is 1000)
%     plot        true or false (default is true)
%     ploteach    plot each time step (default is 10)
%     plotchi     {false} or true. If true, the plot will show the
%                 chi-transformed profile.
%     gifname     By default empty, but if filename is provided, a gif file
%                 will be written to the disk. 
%     gifopts     Structure array with gif options for the function gif
%                 (see help gif)
%     ylim        Two-element vector with minimum and maximum of y-axis if
%                 'plot' is true
%
% Example
%
%     % Run with standard values
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     A   = flowacc(FD);
%     S   = klargestconncomps(STREAMobj(FD,'minarea',1000));
%     z   = fastscape(S,DEM,A);    
%
%     % Variable boundary conditions defined as rates
%
%     bct = table;
%     bct.t = [0 10000 50000 800000]';
%     bct.dz = [0 -0.01 0 0]';
%     z   = fastscape(S,DEM,A,'k',1e-5,'plot',true,...
%                     'bc',bct,'bctype','rate','tspan',800000,...
%                     'u',0.003,'ploteach',20,'dt',500);
%
%     % Variable boundary conditions defined as elevations. Here, we 
%     % simulate a sudden drop in elevation between 200 and 200.5 ky by 
%     % 200 m.
%
%     bct = table;
%     bct.t = [0 200000 200500 400000]';
%     zb  = DEM.Z(streampoi(S,'outlet','ix'));
%     bct.z = [zb zb zb-200 zb-200]';
%     z   = fastscape(S,DEM,A,'k',1e-5,'plot',true,...
%                     'bc',bct,'bctype','elev','tspan',400000,...
%                     'u',0.003,'ploteach',5,'dt',500);
%
%     % And here are multiple drops in elevation
%     bct = table;
%     bct.t = (0:500000)';
%     bct.z = zeros(size(bct.t));
%     bct.z(mod(bct.t,50000)==0) = -50;
%     bct.z = zb + cumsum(bct.z);
%     z   = fastscape(S,DEM,A,'k',1e-5,'plot',false,...
%                     'bc',bct,'bctype','elev','tspan',500000,...
%                     'u',0.003,'dt',50);
%
%     % You can also set a function to control boundary conditions. Here,
%     % we set boundaries as rates of baselevel drop 
%     bcfun = @(t) (sin(t/1e4) - 1)/1000 * 2;
%     z   = fastscape(S,DEM,A,'k',1e-5,'plot',true,...
%                     'bc',bcfun,'bctype','rate','tspan',500000,...
%                     'u',0.003,'dt',50);
%
%     % Here's another example that shows how a knickpoint consumes another
%     % if n ~= 1.
%     ST  = trunk(S);
%     bct = table;
%     bct.t = [0 200000 200500 220000 220500 5e5]';
%     zb  = DEM.Z(streampoi(S,'outlet','ix'));
%     bct.z = [zb zb zb-50 zb-50 zb-200 zb-200]';
%     z   = fastscape(ST,DEM,A,'k',4e-5,'plot',true,...
%                     'bc',bct,'bctype','elev','tspan',5e5,...
%                     'u',0.003,'dt',50, 'n',1.3,'ploteach',50);
%
%
%
% References
%
%     Braun, J. and Willett, S. D.: A very efficient O(n), implicit and
%     parallel method to solve the stream power equation governing fluvial
%     incision and landscape evolution, Geomorphology, 180–181, 170–179,
%     doi:10.1016/j.geomorph.2012.10.008, 2013.
% 
%     Hergarten, S., Neugebauer, H.J., 2001. Self-organized critical 
%     drainage networks. Phys.Rev. Lett. 86, 2689–2692
%
%     Hergarten, S., 2002. Self organised criticality in Earth Systems. 
%     Springer, Heidelberg.
%
%     Campforts, B., Schwanghart, W., Govers, G. (2017): Accurate 
%     simulation of transient landscape evolution by eliminating numerical 
%     diffusion: the TTLEM 1.0 model. Earth Surface Dynamics, 5, 47-66. 
%     [DOI: 10.5194/esurf-5-47-2017]
%
% See also: STREAMobj, gif
%
% Authors: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de) and
%          Benjamin Campforts.
% Date: 12. March, 2024


gifopts.DelayTime = 1/15;
gifopts.LoopCount = inf;
% gifopts.frame     = gcf; % commented out because it will open a window
% even in the case if no plot is wanted.
gifopts.overwrite = false;

%% Input parsing
p = inputParser;
addRequired(p,'S',@(x) isa(x,'STREAMobj'))
addRequired(p,'z')
addRequired(p,'a')

addParameter(p,'convertpx',true)
addParameter(p,'uplift',0.001) % 1 mm y^-1 default
addParameter(p,'k',1e-5)
addParameter(p,'m',0.5)
addParameter(p,'n',1)
addParameter(p,'tspan',100000)
addParameter(p,'dt',1000)
addParameter(p,'bc',[])
addParameter(p,'bctype','rate')
addParameter(p,'plot',true)
addParameter(p,'ploteach',10)
addParameter(p,'plotchi',0)
addParameter(p,'gifname','')
addParameter(p,'gifopts',gifopts)
addParameter(p,'ylim',[],@(x) isempty(x) || ((numel(x) == 2) && (x(1)< x(2))))
parse(p,S,z,a,varargin{:})

% STREAMobj
S = p.Results.S;
% Elevations
z = ezgetnal(S,p.Results.z);
z = imposemin(S,z);

% Set y limits if plot is required
if isempty(p.Results.ylim) 
    ylimauto = true;
else
    ylimauto = false;
    yl = p.Results.ylim;
end

% Upstream areas
a = ezgetnal(S,p.Results.a);
if p.Results.convertpx
    a = a*S.cellsize^2;
end

% Uplift
u = ezgetnal(S,p.Results.uplift);
% Erodibility K
k = ezgetnal(S,p.Results.k);
% Stream power parameter
m = p.Results.m;
n = p.Results.n;
% Simulation time and timestep
tspan = p.Results.tspan;
dte   = p.Results.dt;
% Plot?
plotit   = p.Results.plot;
ploteach = p.Results.ploteach;
plotchi  = p.Results.plotchi;
% Write gif
writegif = p.Results.gifname;
gifopts  = p.Results.gifopts;


%FASTSCAPE1D 1D implementation of Braun and Willett 2013 implicit scheme
ix  = S.ix;
ixc = S.ixc;
d   = S.distance;
dx_ixcix = d(ix)-d(ixc);

% Calculate K*A^m
ar     = a;
a      = k.*(a.^m);

% calculate timesteps
dte    = 0:dte:tspan;
if dte(end) <= tspan
    dte(end+1) = tspan;
end
dte    = diff(dte);

% get timeseries of boundary conditions
outlet = streampoi(S,'outlet','logical');
outletix = find(outlet);
bctype = validatestring(p.Results.bctype,{'rate','elev'});
if isempty(p.Results.bc)
    zb = repmat(z(outletix),1,numel(dte));
elseif isnumeric(p.Results.bc) && strcmp(bctype,'rate')
    zb = z(outletix) + p.Results.bc*cumsum(dte);
elseif isnumeric(p.Results.bc) && strcmp(bctype,'elev')
    zb = interp1([0;tspan],...
                 [z(outletix)'; repmat(p.Results.bc,1,numel(outletix))],...
                 cumsum(dte)');
    zb = zb';
elseif istable(p.Results.bc)
    BCT = p.Results.bc;
    if BCT.t(end) < tspan
        BCT.t(end) = tspan;
        warning('TopoToolbox:fastscape',...
            ['The last element in the time column of the boundary\n' ...
             'table was set to ' num2str(tspan)]);
    end

    dtc = cumsum([0 dte]);
    switch lower(bctype)       
        case 'elev'                     
            zb = interp1(BCT.t,BCT.z,dtc(:));
            zb = zb';
        case 'rate'
            zc = cumtrapz(BCT.t,BCT.dz);
            zb = z(outletix) + interp1(BCT.t,zc,dtc(:));
            zb = zb';      
    end
    
elseif isa(p.Results.bc,'function_handle')
    if strcmp(bctype,'rate')
        dtc = cumsum([0 dte]);
        zb = z(outletix) + cumtrapz(dtc,p.Results.bc(dtc));
    else
        zb = z(outletix)*0 + cumsum(p.Results.bc(dte));
    end
else
    error('Cannot handle boundary conditions')
end

if plotit
    if plotchi == 0
        d  = S.distance;
    elseif plotchi == 1
        d  = chitransform(S,ar,"mn",m/n);
    end
    hh = plotdz(S,z,'color','r','distance',d); 
    if ~ylimauto
        ylim(yl)
    end
    hold on
    if plotchi
        xlabel('\chi [m]')
    end

    if ~isempty(writegif)
        gif(writegif,gifopts)
    end
end

t = 0;
plotcounter = 0;

% linear stream power model (n==1)
if n == 1
    for titer = 1:numel(dte)
        t = t+dte(titer);
        % adjust baselevel
        z(outlet) = zb(:,titer);
        % z(outlet) = z(outlet)-baseleveldrop*dte(titer)/tspan;
        % add uplift
        z(~outlet) = z(~outlet) + u(~outlet)*dte(titer);
        
        for r = numel(ix):-1:1
            tt       = a(ix(r))*dte(titer)/(dx_ixcix(r));
            z(ix(r)) = (z(ix(r)) + z(ixc(r))*tt)./(1+tt);
        end
        
        if plotit && (mod(titer,ploteach)==0 || titer==numel(dte))
            plotcounter = plotcounter + 1;
            if plotcounter > 10
                plotcounter = 10;
                delete(h(1))
                h(1:9) = h(2:10);
            end            
            h(plotcounter) = plotdz(S,z,'color','k','distance',d);
            if plotchi
                xlabel('\chi [m]')
            end
            if plotcounter > 2
                for rr = 1:plotcounter
                    h(rr).Color = [0 0 0 rr/(plotcounter*2)];
                end
            end
            if plotcounter == 1
                ht = text(0.05,0.9,['t = ' num2str(t) ' yrs'],'Units','normalized');
            else
                ht.String = ['t = ' num2str(t) ' yrs'];
            end

            drawnow
            if ~isempty(writegif)
                gif
            end

        end
    end
    
elseif n ~= 1
    
    for titer = 1:numel(dte)
        t = t+dte(titer);
        z(outlet) = zb(:,titer);
        % adjust baselevel
        % z(outlet) = z(outlet)-baseleveldrop*dte(titer)/tspan;
        % add uplift
        z(~outlet) = z(~outlet) + u(~outlet)*dte(titer);
        
        for r = numel(ix):-1:1
            dx      = dx_ixcix(r);
            tt      = a(ix(r))*dte(titer)/dx;
            % z_t
            zt      = z(ix(r));
            % z_(t+dt) of downstream neighbor
            ztp1d   = z(ixc(r));
            % dx
 
            % initial value for finding root
            if ztp1d < zt
                ztp1    = newtonraphson(zt,ztp1d,dx,tt,n);
            else
                ztp1    = zt;
            end
            
            if ~isreal(ztp1) || isnan(ztp1)
                disp('Non real solutions converted to real')
                ztp1=real(ztp1);
            end
            z(ix(r))=ztp1;
            
        end
        if plotit && (mod(titer,ploteach)==0 || titer==numel(dte))
            plotcounter = plotcounter + 1;
            if plotcounter > 10
                plotcounter = 10;
                delete(h(1))
                h(1:9) = h(2:10);
            end            
            h(plotcounter) = plotdz(S,z,'color','k','distance',d);
            if plotchi
                xlabel('\chi [m]')
            end


            if plotcounter > 2
                for rr = 1:plotcounter
                    h(rr).Color = [0 0 0 rr/(plotcounter*2)];
                end
            end
            if plotcounter == 1
                ht = text(0.05,0.9,['t = ' num2str(t) ' yrs'],'Units','normalized');
            else
                ht.String = ['t = ' num2str(t) ' yrs'];
            end

            drawnow
            if ~isempty(writegif)
                gif
            end
        end
    end
end

    function ztp1 = newtonraphson(zt,ztp1d,dx,tt,n)
        
        tempz   = zt;
        tol = inf;
        
        while tol > 1e-3
            % iteratively approximated value
            tempdz = tempz-ztp1d;
            ztp1  =  tempz - (tempz-zt + ...
                (tt*dx) * ...
                (tempdz./dx)^n) / ...
                (1+n*tt*(tempdz./dx)^(n-1));
            tol   = abs(ztp1-tempz);
            tempz = ztp1;
        end
    end

if plotit
    hold off
end

end

%% Easy handling of nals.
function z = ezgetnal(S,z)
if isa(z,'GRIDobj')
    z = double(getnal(S,z));
elseif isnal(S,z)
    z = double(z);
elseif isnumeric(z) && isscalar(z)
    z = getnal(S) + double(z);
else
    error('Cannot handle input')
end
end

