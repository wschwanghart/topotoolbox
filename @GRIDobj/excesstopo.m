function DEM = excesstopo(DEM,maxgradient,varargin)

%EXCESSTOPO Calculate excesstopography using a new algorithm
%
% Syntax
%
%     DEM2 = excesstopo(DEM,maxgradient)
%
% Description
%
%     This function calculates excesstopography in a different way then the
%     function GRIDobj/excesstopography. It enables spatially variable
%     maximum gradients.
%
% Input arguments
%
%     DEM           GRIDobj
%     maxgradient   maximum allowed gradient (scalar, matrix same size as
%                   DEM.Z or GRIDobj)
%
% Output arguments
%
%     DEM2          Reconstructed DEM
%
%                  
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 9. February, 2023

if nargin == 2
    nmax = 1000;
else
    nmax = varargin{1};
end

if isnumeric(maxgradient) && isscalar(maxgradient)
    DEM.Z = excesstopoUniSlope(DEM.Z,DEM.cellsize,maxgradient,nmax);
else
    if isa(maxgradient,'GRIDobj')
        validatealignment(maxgradient,DEM)
        maxgradient = maxgradient.Z;
    else
        assert(isequal(DEM.size,size(maxgradient)))
    end
    DEM.Z = excesstopoVarSlope(DEM.Z,DEM.cellsize,maxgradient,nmax);
end