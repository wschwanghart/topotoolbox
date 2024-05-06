function [Rb,Ns,s] = bifurcationratio(S)

%BIFURCATIONRATIO Calculate the bifurcation ratio of a STREAMobj
%
% Syntax
%
%     [Rb,Ns,s] = bifurcationratio(S)
%
% Description
%
%     The bifurcation ratio is the average ratio between the number of
%     streams of a given order and the number of streams in the next higher
%     order. This function calculates the bifurcation ratio for a stream
%     network stored in a STREAMobj S. The calculation is not performed
%     for individual basins, but for the entire network. If you want to
%     calculate the Rb for each basin, please see the example below.
%
% Input arguments
%
%     S      STREAMobj
%     
% Output arguments
%
%     Rb     bifurcation ratio (scalar)
%     Ns     Number of streams of a given order, so that Ns(1) is the
%            number of first order streams, Ns(2) is the number of second
%            order streams, ...
%     s      node-attribute list with streamorders
%
% Example: Calculate Rb for several drainage basins
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',500);
%     S = klargestconncomps(S,3);
%     CS = STREAMobj2cell(S);
%     Rb = cellfun(@(S) bifurcationratio(S),CS);
%
%     % Plot the bifurcation ratio as a node-attribute list
%     Rbnal = getnal(S);
%     for r = 1:numel(CS) 
%         Rbnal = nal2nal(S,CS{r},getnal(CS{r})+Rb(r),Rbnal); 
%     end
%     plotc(S,Rbnal)
%
% See also: STREAMobj, STREAMobj/streamorder
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 6. May, 2024

% Calculate streamorder
s = streamorder(S);

% Maximum streamorder
maxs = max(s);

% If there are only first-order streams, Rb is not defined.
if maxs == 1
    Rb = nan;
    [~,Ns] = conncomps(S);
    return
end

% Preallocate array for number of streams of a given order
Ns   = zeros(maxs,1);

for r = 1:maxs
    % Extract streams with order r
    S2 = subgraph(S,s == r);
    % Calculate the number of connected components 
    [~,Ns(r)] = conncomps(S2);
end

% Calculate Rb
Rb = 1/(maxs-1)*sum(Ns(1:end-1)./Ns(2:end));
