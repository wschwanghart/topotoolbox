function varargout = cdftest(P,varargin)

%CDFTEST Kolmogorov-Smirnov test based on the CDF of a covariate
%
% Syntax
%
%     t = cdftest(P)
%     t = cdftest(P,'pn',pv)
%     [t,p,ks2stat] = ...
%
% Description
%
%     cdftest performs a test of CSR using the Kolmogorov-Smirnov
%     statistics (KS-Test).
%
%     The test returns t = 0, if the null hypothesis cannot be rejected
%     that the point pattern in P is completely random.
%
% Input arguments
%
%     P      PPS
%     
%     Parameter name/value pairs
% 
%     'covariate'    P.S.distance
%     other parameters see kstest2
%
% Output arguments
%
%     see kstest2
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S,1);
%     P = PPS(S,'runif',100,'z',DEM);
%     t = cdftest(P); % should return 0
%
% See also: PPS/ecdf, PPS/quadratcount, PPS/rhohat
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 21. November, 2022

p = inputParser;
p.FunctionName = 'PPS/cdftest';

addParameter(p,'covariate',P.S.distance,@(x) isnal(P.S,x) | isa(x,'GRIDobj'));
addParameter(p,'alpha',0.05)
addParameter(p,'tail','unequal')

parse(p, varargin{:})

x1 = getcovariate(P,p.Results.covariate);
x2 = x1(P.PP);

varargout = cell(1,max(nargout,1));

[varargout{:}] = kstest2(x1,x2,'Alpha',p.Results.alpha,'Tail',p.Results.tail);