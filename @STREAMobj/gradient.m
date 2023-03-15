function s = gradient(S,DEM,varargin)

%GRADIENT along-stream gradient
%
% Syntax
%
%     s = gradient(S,DEM)
%     s = gradient(S,DEM,pn,pv,...)
%     s = gradient(S,z)
%     s = gradient(S,z,pn,pv,...
%
% Description
%
%     GRADIENT calculates the stream slope for each node in the stream
%     network S based on the associated digital elevation model DEM.
%     Different methods (see parameter name value pairs) can be applied to
%     calculate the slope. The function returns a node-attribute list.
%
% Input arguments
%
%     S    instance of STREAMobj
%     DEM  digital elevation model (class: GRIDobj)
%     z    node attribute list
%
%     parameter name/value pairs {default}
%
%     'unit' string
%     {'tangent'} 'degree' 'radian' 'percent' 'sine'
%
%     'method' string
%     {'forward'} is the same as steepest downward gradient. 
%     'robust' calculates the gradient based on a minimum, vertical drop
%     (see 'drop'-option). 'robust' ensures that all gradients are nonzero,
%     however gradient may not be computed close to the edges.
%     
%     'drop' scalar {10}
%     only applicable when 'method' is set to 'robust'. Determines the
%     minimum vertical drop that a stream must have to calculate the
%     gradient.
%
%     'imposemin' {false} or true
%     minima imposition to avoid negative slopes (see imposemin)
%
% Output arguments
%
%     s      stream gradient as node-attribute list
%
% 
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     S  = STREAMobj(FD,A>1000);
%     S  = klargestconncomps(trunk(S));
%     g  = gradient(S,DEM,'method','robust');
%     subplot(2,1,1); plotdz(S,DEM);
%     subplot(2,1,2); plotdz(S,g); ylabel('Gradient [-]')
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 3. May, 2013

narginchk(2,inf)

% Parse Inputs
p = inputParser;         
p.FunctionName = 'gradient';

validunits = {'tangent' 'degree' 'radian' 'percent' 'sine'};
validmethods = {'forward' 'centered' 'robust'};

addParamValue(p,'unit','tangent',@(x) ischar(validatestring(x,validunits)));
addParamValue(p,'method','forward',@(x) ischar(validatestring(x,validmethods)));
addParamValue(p,'drop',10,@(x) isscalar(x) && x>0);
addParamValue(p,'imposemin',false,@(x) isscalar(x));

parse(p,varargin{:});

% get node attribute list with elevation values
if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    z = double(getnal(S,DEM));
elseif isnal(S,DEM);
    z = double(DEM);
else
    error('Imcompatible format of second input argument')
end

% if imposemin
if p.Results.imposemin
    for r = 1:numel(S.ix);
        z(S.ixc(r)) = min(z(S.ix(r)),z(S.ixc(r)));
    end
end

d = hypot(S.x(S.ix)-S.x(S.ixc),S.y(S.ix)-S.y(S.ixc));

switch validatestring(p.Results.method,validmethods);
    case 'forward'
        s       = zeros(size(S.x));
        s(S.ix) = (z(S.ix)-z(S.ixc))./d;
    case 'centered'
        error('centered method is not implemented yet')
    case 'robust'
        
        % fast downstream indexing
        ixcix  = zeros(numel(S.x),1);
        ixcix(S.ix) = 1:numel(S.ix);
        
        % preallocate slope
        s = nan(size(S.x));
        drop = p.Results.drop;
        
        for r = 1:numel(S.ix);
            
            zz = z(S.ix(r));
            dd = 0;
            
            c  = r;
            
            % first go until a cell height difference is less than the drop
            while (c ~= 0) && (zz-z(S.ixc(c)) < drop);
                dd = dd+d(c);
                c = ixcix(S.ixc(c));
            end
            
            if c == 0;
                s(S.ix(r)) = nan;
            else
%                 z1 = z(S.ix(c));
                z2 = z(S.ixc(c));
                
                s(S.ix(r)) = (zz-z2)./(d(c) + dd);
                % lower the subsequent cell a little (check, if this is ok)
                z(S.ixc(r)) = min(z(S.ix(r)) - s(S.ix(r))*d(r),z(S.ixc(r)));
%                 
%                 z1  = -(zz-z(S.ix(c))-drop);
%                 z2  = zz-z(S.ixc(c))-drop;
%                 
%                 ddd = z2/z1 *d(c) + dd;
%                 
%                 s(S.ix(r)) = drop./ddd;
            end
        end
end
                
                
                
                
            
            



switch validatestring(p.Results.unit,validunits);
    case 'tangent'
        % do nothing
    case 'degree'
        s = atand(s);
    case 'radian'
        s = atan(s);
    case 'sine'
        s = sin(atan(s));
    case 'percent'
        s = s*100;
end
