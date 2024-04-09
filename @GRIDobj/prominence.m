function [x,y,p] = prominence(DEM,tol,minelev)

%PROMINENCE Calculate the prominence of mountain peaks
%
% Syntax
%
%     [x,y,p] = prominence(DEM,tol)
%
% Description
%
%     This function calculates the prominence of peaks in a DEM. The
%     prominence is the minimal amount one would need to descend from a
%     peak before being able to ascend to a higher peak.
%
%     The function uses image reconstruct (see function imreconstruct) to
%     calculate the prominence. It may take a while to run for large DEMs.
%     The algorithm iteratively finds the next higher prominence and stops
%     if the prominence is less than the tolerance, the second input
%     parameter to the function.
%
%     The function opens a waitbar which let's you quit the operation. The
%     list of coordinates and prominence values derived until then, are
%     returned as output arguments.
%
% Input arguments
%
%     DEM    Digital elevation model (GRIDobj)
%     tol    tolerance in [m] with the minimum prominence.
%
% Output arguments
%
%     x,y    coordinate values of peaks
%     p      prominence (the first element in the vector contains the
%            highest peak for which prominence (given above definition)
%            cannot be calulated.
%
% Example
%
%     DEM = readexample('taiwan');
%     [x,y,p] = prominence(DEM,200);
%     imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false);
%     hold on
%     h = bubblechart(x,y,p);
%     bubblelim([200 2000])
%     bubblesize([2 20])
%     h.MarkerEdgeColor = 'k';
%     bubblelegend('location','southeast')
%
% See also: GRIDobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 29. June, 2021

narginchk(2,3)

if any(isnan(DEM))
    DEM.Z(isnan(DEM.Z)) = 0;
end

if nargin < 3
    minelev = -inf;
end

P = GRIDobj(DEM)+min(DEM);
P2 = P;
C  = P;

counter = 1;
p = inf;

f = waitbar(0,'1','Name','Calculating prominence...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);

elev = inf;

while p(end) > tol 
     
    [p(counter),ix(counter)] = max(DEM-P);
    
    waitbar(min(max(tol/p(end),0),1),f,sprintf('%5.2f',p(end)))
    
    P.Z(ix) = DEM.Z(ix);
    P.Z = imreconstruct(P.Z,DEM.Z);
    counter = counter + 1;
    elev = DEM.Z(ix);
    
    % Check for clicked Cancel button
    if getappdata(f,'canceling')
        break
    end
    
end

delete(f);
[x,y] = ind2coord(DEM,ix);
p = p(:);

