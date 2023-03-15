function [ixchannel,d,x,y] = flowpathextract(FD,ixchannel,A,stopcrit)

%FLOWPATHEXTRACT extract linear indices of a single flowpath in a DEM
%
% Syntax
%
%     [IX,distance] = flowpathextract(FD,channelstart)
%     [IX,distance] = flowpathextract(FD,channelstart,A)
%     [IX,distance] = flowpathextract(FD,channelstart,A,stopcrit)
%     [IX,distance,x,y] = ...
%
% Description
%
%     flowpathextract finds the linear indices of a flow path starting at
%     channelstart based on flow directions given an instance of FLOWobj.
%     If called repeatedly, consider setting the fastindexing property of
%     the FLOWobj to true (set.fastindexing(FD,true)). flowpathextract will
%     be much faster than.
%
%     If called with three or four input arguments, the flowpath will be 
%     derived in upstream direction. In contrast to downstream flowpath
%     extraction, usually a large number of possible flowpaths can be 
%     derived in upstream direction. The GRIDobj A provides weights according
%     to which the path is extracted. A is usually a flow accumulation grid
%     obtained from the function FLOWobj/flowacc and thus, the flowpath will
%     follow the main river trunk to the watershed. stopcrit allows to set
%     a stop criterion at which flow path extraction terminates and will.
%     usually be a minimum upslope area in the same units as A. 
%     
%     Be advised that the algorithm uses rand to avoid the extraction of 
%     multiple paths and thus flow paths may vary for same inputs. The 
%     intercell variability in the weight raster should thus be at best
%     larger than one. Otherwise, the flow path will be highly random.
%
% Input
%
%     FD            flow direction object (FlowDirObj)
%     channelstart  linear index (scalar) where channel starts
%     A             weight grid (GRIDobj) as obtained by flowacc. If
%                   supplied, the flow path will be extracted in upstream
%                   direction
%     stopcrit      value in A at which upstream flow path extraction 
%                   terminates (e.g. minimum flow accumulation, same units
%                   as in A)
%
% Output
%
%     IX            linear index of channel pixels ordered in downstream
%                   order
%     distance      horizontal distance along flowpath
%     x,y           coordinate vectors of flowpath
%
% Example
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     [IXc,d] = flowpathextract(FD,473962);
%     plot(d,DEM.Z(IXc))
%
% See also: FlowDirObj, GRIDobj/coord2ind
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 26. January, 2013


narginchk(2,4)
validateattributes(ixchannel,{'numeric'},{'scalar','integer'},'flowpathextract','channelstart',2)

if nargin == 2;
    type = 'downstream';
else
    type = 'upstream';
    if nargin == 3
        stopcrit = -inf;
    end
end
                
switch type
    case 'downstream'
        
        
        if ~FD.fastindexing;
            FD.fastindexing = true;
        end
        
        r = 1;
        
        ixcix = FD.ixcix;
        ixc   = FD.ixc;
        
        while ixcix(ixchannel) ~= 0
            r = r+1;
            ixchannel(r) = ixc(ixcix(ixchannel(end)));
        end
        
    case 'upstream'
        
        % evaluate the costs to move from one pixel to neighboring,
        % upstream pixels. This is a little tricky since each pixel can
        % have up to 8 upstream neighbors
        
        % value of upstream neighbor
        a = A.Z(FD.ix);
        I = a>stopcrit;
        
        % pruning the network using the threshold criterion reduces the
        % number of links
        ix = FD.ix(I);
        ixc = FD.ixc(I);
        a  = double(a(I)) + rand(numel(ix),1);
        
        % minimal example
        %
        % ix = [1       ixc = [2         a = [.5
        %       2              4              1 
        %       3              4              3
        %       4]             5]             6]
        %
        % should yield
        %
        %  ix  ixc
        %   1    2 
        %   3    4
        %   4    5
        %
        % multiple edges from a have been reduced to one edge
        % according to the maximum value in a
        %
        
        Amax = zeros(FD.size);
        for r = 1:numel(ix)
             Amax(ixc(r)) = max(Amax(ixc(r)),a(r));
        end
        
        I = Amax(ixc) == a;
        
        ix = ix(I);
        ixc = ixc(I);
        
        % fast indexing in upstream direction
        ixixc = zeros(FD.size,'uint32');
        ixixc(ixc) = uint32(1):uint32(numel(ixc));
        
        r = 1;
        while ixixc(ixchannel) ~= 0
            r = r+1;
            ixchannel(r) = ix(ixixc(ixchannel(end)));
        end
        
        ixchannel = ixchannel(end:-1:1);
        
        
end

ixchannel = ixchannel(:);

if nargout > 1;
    d = getdistance(ixchannel(1:end-1),ixchannel(2:end),FD.size,FD.cellsize);
    d = [0; cumsum(d)];
end

if nargout > 2;
    C = copy2GRIDobj(FD);
    [x,y] = ind2coord(C,ixchannel);
end
