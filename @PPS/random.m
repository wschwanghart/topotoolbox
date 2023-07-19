function Psim = random(P,mdl,covariates)

%RANDOM random realisation of a (inhomogeneous) point process
%
% Syntax
%
%     Psim = random(P)
%     Psim = random(P,mdl)
%     Psim = random(P,mdl,covariates)
%
% Description
%
%     random creates a random realisation of a point pattern based on the 
%     loglinear model in mdl. 
%
%     random(P) returns a point pattern with homogeneous intensity derived
%     from P.
%
%     random(P,mdl) returns a point pattern with inhomogeneous intensity
%     derived from the GeneralizedLinearModel mdl derived from
%     fitloglinear.
%
%     random(P,mdl,covariates) returns a point pattern with inhomogeneous
%     intensity based on mdl and covariates.
%
% Input arguments
%
%     P           instance of PPS
%     mdl         loglinear model
%     covariates  node-attribute lists of independent variables in the 
%                 model. The lists must be compatible with mdl.
%
% Output arguments
%
%     Psim        new instance of PPS with random points
%
% 
% See also: PPS, PPS/fitloglinear, PPS/simulate
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019


if nargin == 1
    Psim = PPS(P.S,'rpois',intensity(P));
    P.PP = Psim.PP;
    Psim = P;
    return
end

if nargin == 2
    covariates = mdl.Variables(:,1:end-1);
end

ysim = random(mdl,covariates,'offset',mdl.Offset);
Psim = P;
Psim.PP = find(ysim);
