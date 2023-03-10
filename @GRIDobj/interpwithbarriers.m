function ZI = interpwithbarriers(DEM,x,y,z,varargin)

%INTERPWITHBARRIERS Laplace interplation with barriers
%
% Syntax
%
%     ZI = interpwithbarriers(DEM,x,y,z)
%     ZI = interpwithbarriers(DEM,x,y,z,'pn',pv)
%
% Description
%
%     interpwithbarriers implements a laplacian interpolation with
%     barriers. While laplacian interpolation generates a smoothly varying
%     surface (the ridge parameter) enables increasing the smoothness,
%     there may be sharp local differences, i.e. along fault systems or 
%     boundaries. Note that together with a ridge parameter > 0, the
%     interpolator is not exact, i.e., the surface passes not through the
%     observations.
%
% Input arguments
%
%     DEM     GRIDobj that is used to provide a geometry
%     x,y     Coordinate vectors of observations
%     z       Observations
%     
%     Parameter name/value pairs
%
%     'barriers'   structure array of lines (mapping structure) or logical 
%                  GRIDobj with ones indicating lines (see line2GRIDobj).
%     'ridge'      0 by default. Larger values lead to higher smoothing
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%      
%     [x,y] = randomsample(DEM,100);
%     
%     % create a barrier line (cross)
%     [xx,yy] = getcoordinates(DEM);
%     BARRIER = GRIDobj(DEM);
%     
%     BARRIER.Z(150:450,600) = true;
%     BARRIER.Z(300,300:900) = true;
%     
%     z = rand(size(x));
%     z(x > xx(600)) = z(x > xx(600)) + 1;
%     DEM = GRIDobj(DEM)+3;
%     % DEM.Z(2:end-1,2:end-1) = nan;
%     
%     Z = interpwithbarriers(DEM,x,y,z,'barrier',BARRIER,'ridge',20);
%     
%     figure
%     surf(Z)
%     exaggerate(gca,700)
%     camlight
%     hold on
%     stem3(x,y,z,'k')
%     hold off
%
% See also: GRIDobj/resample, GRIDobj/reclabel, accumarray
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 10. March, 2023

p = inputParser;
p.FunctionName = 'GRIDobj/laplaceinterp';
addParameter(p,'barriers',[]);
addParameter(p,'ridge',0)
addParameter(p,'lb',[])
addParameter(p,'ub',[])
parse(p,varargin{:});
    
% Ensure column vectors
x = x(:);
y = y(:);
z = z(:);

% Make sure that they have the same size
assert(isequal(size(x),size(y)) & isequal(size(x), size(z)), ...
    'x, y and z must have equal size.' )

ZI = aggregate(DEM,[x y z],@mean);

G = ~GRIDobj(DEM,'logical');
[ic,icd] = ixneighbors(G.Z,G.Z,4);

if ~isempty(p.Results.barriers)
    
    if isstruct(p.Results.barriers)
        F = line2GRIDobj(DEM,p.Results.barriers);
    else
        F = p.Results.barriers;
        validatealignment(F,DEM)
    end

    I = F.Z(ic) | F.Z(icd);
    ic(I) = [];
    icd(I) = [];
end

n = numel(DEM.Z);
A = sparse(ic,icd,1,n,n);

inan = isnan(ZI.Z);

s  = 1./(sum(A,2));
s(~inan(:)) = 0;
A = spdiags(s,0,n,n)*A;

% fill main diagonal with ones
A = speye(n)-A;

% solve
ZI.Z(inan) = 0;
if p.Results.ridge == 0

    ZI.Z = reshape(A\ZI.Z(:),DEM.size);

else

    Asd  = sparse(ic,icd,1,n,n);
    Asd  = spdiags(sum(Asd,2),0,n,n) - Asd;
    Asd  = p.Results.ridge*Asd;
    
    if isempty(p.Results.lb) && isempty(p.Results.ub) 

        zi   = [A;Asd]\[ZI.Z(:);zeros(n,1)];
    
    else

        C = [A;Asd];
        d = [ZI.Z(:);zeros(n,1)];
        
        b = p.Results.lb.Z(:);
        zi = -lsqlin(C,-d,speye(n),-b);

    end
    ZI.Z = reshape(zi,DEM.size);
    
end

% Deal with the barrier which is nan or zero
ZI.Z(F.Z>0) = nan;
ZI = inpaintnans(ZI);
    