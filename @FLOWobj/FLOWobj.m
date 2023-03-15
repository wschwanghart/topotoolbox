classdef FLOWobj
    
%FLOWOBJ create flow direction object
%
% Syntax
%
%     FD = FLOWobj(DEM)
%     FD = FLOWobj(DEM,'pn','pv',...)
%     FD = FLOWobj(DEM,'multi')
%     FD = FLOWobj(DEM,'Dinf')
%     FD = FLOWobj(M,'cellsize',cs,'size',siz)
%
% Description
%
%     FLOWobj creates a flow direction object which can be used to
%     calculate terrain attributes such as flow accumulation, drainage
%     basin delineation, flow path extraction, etc. In contrast to the M
%     matrix returned by the function flowdir, the flow direction object
%     requires much less memory and the same routines are usually evaluated
%     much faster.
%
%     FD can be calculated from the original digital elevation model or
%     converted from an existing flow direction matrix. FD contains a
%     topological ordering of the DEM nodes.
%
%     FLOWobj(DEM,'multi') derives multiple flow direction (MFD) and
%     FLOWobj(DEM,'dinf') derives D infinity according to Tarboton's (1997)
%     method (Eddins 2016). Note that this option precludes any other 
%     options to be set. Thus, this method needs preprocessing, e.g. using 
%     fillsinks(DEM).
%
% Input arguments
%
%     DEM    digital elevation model (Class: GRIDobj)
%     M      flow direction matrix as returned by flowdir, flowdir_single
%            (function of TopoToolbox 1)
%     
% Parameter name/value pairs   {default}
%
%  Applicable only, if calculated from GRIDobj
%
%     'preprocess' --  {'carve'}, 'fill', 'none'
%            set DEM preprocessing that determines flow behavior in
%            topographic depressions and flat areas 
%     'sinks' -- logical matrix same size as dem
%            true values in sinks are treated as sinks in the digital
%            elevation model and are not filled or carved, if the
%            preprocessing option fill or carve are chosen.
%     'internaldrainage' -- {false} or true
%            set this parameter value to true if flow directions should be
%            derived in the lowest, flat regions of internal drainage 
%            basins. By default, this parameter is set to false since this
%            information is usually not required and flow paths will stop
%            when entering flat, internally drained sections.
%     'cweight' --  scalar {1}
%            adjust cost matrix if preprocessing option 'carve' has been
%            chosen. 
%     'verbose' --  {false},true
%            verbose output in the command window to track computational
%            progress. Particularly interesting when working with very
%            large matrices.
%     'mex' --  {'false'},true
%            controls if the mex routines should be called. If true, the mex
%            functions (see function compilemexfiles) must be available on
%            the search path. Using mex functions can increase the speed at
%            which an instance of FLOWobj is constructed.
%
% Applicable only, if calculated from flow direction matrix M
%
%     'cellsize' -- {1}
%            cellsize in x and y direction (scalar double) 
%     'size' --  (1x2 array, e.g. size(dem))
%            the size of the DEM of which the flow direction matrix was 
%            derived from. 
%     'refmat' -- referencing matrix as derived from GRIDobj property refmat
%     'algorithm' -- {'dmperm'} or 'tsort' or 'toposort'
%            determines the algorithm to perform a topological sorting
%            of the vertices of the directed graph in matrix M (this option
%            is left here mainly for evaluation purposes of different 
%            algorithms).
%
% Output
%
%     FD     flow direction object (FLOWobj)
%
% Example 1
%
%     % create and evaluate flow direction object for a DEM
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     A  = flowacc(FD);
%     imageschs(DEM,log(A));
%
% Example 2
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','carve','mex',true);
%     DEM = imposemin(FD,DEM,0.0001);
%     S   = STREAMobj(FD,'minarea',1000);
%     DEMc = DEM;
%     DEMc.Z(S.IXgrid) = DEMc.Z(S.IXgrid)-100; 
%     FD  = FLOWobj(DEMc,'multi');
%     A   = flowacc(FD);
%     imageschs(DEM,log(A));
% 
%
% See also: GRIDobj, STREAMobj
%
% References: 
%
%     Tarboton, D. G. (1997). A new method for the determination of flow 
%     directions and upslope areas in grid digital elevation models. 
%     Water Resources Research, 33(2), 309-319.
%
%     Eddins, S. (2016). Upslope area function. Mathworks File Exchange, 
%     https://www.mathworks.com/matlabcentral/fileexchange/15818-upslope-area-functions
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 02. September, 2017


% Update 
% 2015-06-16: added support for digraph toposort algorithm available with
% R2015b
% 2016-10-07: added support for multiple flow directions
% 2016-11-14: added support for Dinf
% 2017-09-02: default preprocessing option is carve

properties(GetAccess = 'public', SetAccess = 'public')
    size      % size of instance of GRIDobj from which STREAMobj was derived
    type      % flow direction type (single, multi)
    ix        % [edge attribute] topologically sorted nodes (givers)
    ixc       % [edge attribute] topologically sorted nodes (receivers)
    fraction  % [edge attribute] fraction transfered between nodes
    cellsize  % cellsize of the grid (scalar)
    refmat    % 3-by-2 affine transformation matrix (see makerefmat)
    georef    % additional information on spatial referencing
    
end

properties(GetAccess = 'public', SetAccess = 'public')
    fastindexing = false; % set to true to initiate fast indexing
    ixcix     % indexing matrix for fast indexing
end

methods
    function FD = FLOWobj(DEM,varargin)

        if nargin == 0
            % Create an empty FLOWobj
        elseif nargin >= 1 && nargin ~= 2

            % Parse inputs
            p = inputParser;
            p.FunctionName = 'FLOWobj';
            expectedPreProcess = {'none','fill','carve'};
            expectedAlgorithms = {'dmperm', 'toposort','tsort'};

            addRequired(p,'DEM',@(x) issparse(x) || isa(x,'GRIDobj'));
            
            addParamValue(p,'size',[],@(x) isempty(x) || numel(x)==2);
            addParamValue(p,'cellsize',1,@(x) isscalar(x));
            addParamValue(p,'refmat',[]);
            addParamValue(p,'preprocess','carve',@(x) ischar(validatestring(x,expectedPreProcess)));
            addParamValue(p,'tweight',2,@(x) isscalar(x));
            addParamValue(p,'cweight',1,@isnumeric);
            addParamValue(p,'sinks',[],@(x) isa(x,'GRIDobj'));
            addParamValue(p,'verbose',false,@(x) isscalar(x) && islogical(x));
            addParamValue(p,'streams',1);
            addParamValue(p,'weights',[]);
            addParamValue(p,'internaldrainage',false,@(x) isscalar(x) && islogical(x));
            addParamValue(p,'mex',false,@(x) isscalar(x) && islogical(x));
            addParamValue(p,'algorithm','dmperm',@(x) ischar(validatestring(x,expectedAlgorithms)));
            
            addParamValue(p,'type','single');

            parse(p,DEM,varargin{:});

            % required
            DEM        = p.Results.DEM;
            siz        = p.Results.size;
            tweight    = p.Results.tweight;
            preprocess = validatestring(p.Results.preprocess,expectedPreProcess);
            cweight    = p.Results.cweight;
            sinks      = p.Results.sinks;
            verbose    = p.Results.verbose;
            streams    = p.Results.streams;
            weights    = p.Results.weights;
            mexx       = p.Results.mex;
            algo       = validatestring(p.Results.algorithm,expectedAlgorithms);
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if mexx
                % check if compiled mex files are really available
                ext = mexext;
                mexx = exist(['steepestneighbor_mex.' ext],'file')==3 && ...
                       exist(['tsort_mex.' ext],'file')==3;
                if mexx==0
                    warning('TopoToolbox:FLOWobj',...
                            ['cannot find compiled mex functions on the search path. \n'...
                             'FLOWobj continues with using slower m-functions']);
                end
            end
                

            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if issparse(DEM)
                % construct flow direction object from M
                M = DEM;
                FD.cellsize = p.Results.cellsize;
                FD.refmat   = p.Results.refmat;
                
                if ~isempty(siz)
                    FD.size = siz;
                else
                    fprintf(['\n       The size of the DEM was not provided. It will be \n'...
                        '       determined automatically. For large datasets it \n'...
                        '       may be more efficient to provide the size as \n'...
                        '       optional input argument.\n \n']);


                    [ic,icd] = find(M);
                    offsets  = unique(ic-icd);
                    FD.size(1) = max(abs(offsets)-1);
                    FD.size(2) = size(M,1)./FD.size(1); %#ok<CPROP>

                    fprintf('       The size of the DEM was determined to be [%d %d].\n \n',...
                        FD.size(1),FD.size(2));
                end
                
                % topologically sort M and generate list of ordered
                % vertex links

                mtf = any(sum(spones(M),2)>1);

                switch algo
                    case 'tsort'                  
                        if mtf
                            [FD.ix,FD.ixc,FD.fraction] = tsort(M);
                            FD.type = 'multi';
                        else
                            [FD.ix,FD.ixc] = tsort(M);
                            FD.type = 'single';
                        end
                        
                    case {'dmperm','toposort'}
                        switch algo
                            case 'dmperm'
                                [p,~] = dmperm(speye(size(M))-M);
                            case 'toposort'
                                if verLessThan('matlab','8.6')
                                    error('This option requires MATLAB R2015b or later')
                                end
                                G = digraph(M);
                                p = toposort(G);
                        end
                                     
                        if mtf
                            [FD.ix,FD.ixc,FD.fraction] = find(M(p,p));
                            FD.type = 'multi';
                        else
                            [FD.ix,FD.ixc] = find(M(p,p));
                            FD.type = 'single';
                        end
                        p = p(:);
                        FD.ixc = p(FD.ixc);
                        FD.ix  = p(FD.ix);

                end
                
                FD.ix = uint32(FD.ix);
                FD.ixc = uint32(FD.ixc);

                
                
            
            %%  construct flow direction object from DEM
            else
                
                % transfer object properties from DEM to FD
                FD.cellsize = DEM.cellsize;
                FD.refmat   = DEM.refmat;
                FD.georef   = DEM.georef;
                FD.size     = DEM.size;
                FD.type     = 'single';
                nrc         = numel(DEM.Z);
                
                % start preprocessing (e.g. carve or fill)    
                switch lower(preprocess)                   
                    
                    case 'fill'
                        % fill will fill all internally drained basins
                        % (unless specified in the sinks raster). At a
                        % later stage, FLOWobj will find a route along the 
                        % centerline of flat areas.
                        if verbose
                            disp([datestr(clock) ' -- Sink filling'])
                        end
                        
                        if isempty(sinks)
                            DEM = fillsinks(DEM);
                        else
                            DEM = fillsinks(DEM,sinks);
                        end
                        
                    case 'carve'
                        % carve will make use of the topography in
                        % depressions to derive the most realistic flow
                        % paths.
                        if verbose
                            disp([datestr(clock) ' -- Sink filling'])
                        end
                        if isempty(sinks)
                            DEMF = fillsinks(DEM);
                        else
                            DEMF = fillsinks(DEM,sinks);
                        end
                        
                        % By default, weights are calculated as the
                        % difference between the filled and the actual DEM.
                        % There is also a weights option, which is a
                        % GRIDobj with weights. But this is currently
                        % undocumented
                        if isempty(weights)                            
                            D    = DEMF-DEM;
                        else
                            D    = -weights;
                        end
                        D = D.Z;
                        DEM  = DEMF;
                end

                % construct height graph
                % After DEM filling, this code will identify flat sections,
                % sills, and internally drained basins. By default
                % internaldrainage is set to false. That means that FLOWobj
                % will not attempt to route through the lowest region in a
                % internally drained basin (e.g. a flat lake surface).
                if p.Results.internaldrainage
                    [Iobj,SILLSobj,IntBasin] = identifyflats(DEM);
                else
                    [Iobj,SILLSobj] = identifyflats(DEM);
                end
                
                I     = Iobj.Z;
                SILLS = SILLSobj.Z;
                clear Iobj SILLSobj
                
                % calculate sills for internal lake basins. These should be
                % located within the lake to force convergent flows
                if p.Results.internaldrainage
                    % There are various ways to get a lake center pixel
                    % Here we choose the distance transform from outside
                    % the lakes to the inside and take the locations as sills
                    % where the distance is maximum.
                    
                    DD = bwdist(~IntBasin.Z,'e');
                    STATS   = regionprops(IntBasin.Z,DD,'PixelIdxList','PixelValues');
                    if ~isempty(STATS)
                        for r=1:numel(STATS)
                            [~,ixm] = max(STATS(r).PixelValues);
                            STATS(r).MaxIntIX = STATS(r).PixelIdxList(ixm);
                            
                            I(STATS(r).PixelIdxList(1)) = false;
                            SILLS(STATS(r).PixelIdxList(1)) = true;
                        end
                        
                        ixm = [STATS(r).MaxIntIX];
                        I(ixm) = false;
                        SILLS(ixm) = true;
                    end
                    clear InBasin STATS
                    
                    % A slightly faster but less elegant approach
                    % STATS   = regionprops(IntBasin.Z,'PixelIdxList');
                    % for r = 1:numel(STATS);
                    %    I(STATS(r).PixelIdxList(1)) = false;
                    %    SILLS(STATS(r).PixelIdxList(1)) = true;
                    %end
                    % clear IntBasin STATS
                end

                if verbose
                    disp([datestr(clock) ' -- Flat sections identified'])
                end

                % Some more preprocessing required. If the option carve is
                % chosen, we derive here the costs to route through sinks
                switch preprocess
                    case 'carve'
                        CarveMinVal = 0.1;
                        if ~isscalar(cweight)
                            D = (D + cweight);                            
                            D = linscale(D,0,100);
                        end
                        
                        % -- New version -- slightly faster but may require
                        % more memory
%                         STATS = regionprops(I,D,{'PixelIdxList','MaxIntensity','PixelValues'});
%                         PixelValues = cellfun(@(pixval,maxval) (maxval-pixval).^tweight + CarveMinVal,...
%                             {STATS.PixelValues},{STATS.MaxIntensity},...
%                             'UniformOutput',false);
%                         D(vertcat(STATS.PixelIdxList)) = cell2mat(PixelValues);
                        
                        % -- Old version
                        CC = bwconncomp(I);                                
                        for r = 1:CC.NumObjects
                            maxdepth = max(D(CC.PixelIdxList{r}));
                            D(CC.PixelIdxList{r}) = (maxdepth - D(CC.PixelIdxList{r})).^tweight + CarveMinVal;
%                             D(CC.PixelIdxList{r}) = maxdepth - D(CC.PixelIdxList{r});
%                             D(CC.PixelIdxList{r}) = (D(CC.PixelIdxList{r})./maxdepth + CarveMinVal).^tweight;
                        end
                        clear CC

                        % enable that flow in flats follows digitized
                        % streams. Undocumented...
                        if ~isscalar(streams)
                        if islogical(streams)                        
                            D(streams) = CarveMinVal;
                        else
                            D = max(D - D.*min(streams,0.99999),CarveMinVal);
                        end
                        end
                end
               
                % establish the connectivity between sills and flats
                dem = DEM.Z;
                [row,col] = find(SILLS);
                IXsill    = sub2ind(FD.size,row,col);
                rowadd = [-1 -1 0 1 1  1  0 -1];
                coladd = [ 0  1 1 1 0 -1 -1 -1];
                PreSillPixel = [];
                for r = 1:8
                    rowp = row + rowadd(r);
                    colp = col + coladd(r);
                    
                    ValidRowColPair    = rowp>0 & colp>0 & rowp<=FD.size(1) & colp<=FD.size(2);
                    IXPreSill = sub2ind(FD.size,rowp(ValidRowColPair),colp(ValidRowColPair));
                    
                    PreSillPixel = [PreSillPixel;...
                        IXPreSill((dem(IXsill(ValidRowColPair)) == dem(IXPreSill)) & I(IXPreSill))];   %#ok<AGROW>
                    
                end
                
                clear row col IXsill rowadd coladd ValidRowColPair IXPreSill

                % Some more preprocessing if option fill is chosen. Here we
                % derive the costs to route over flats
                I = ~I;
                switch lower(preprocess)
                    case {'fill','none'}
                        D = bwdist(I,'euclidean');
                        mask = inf(FD.size,class(D));
                        mask(I) = 0;
                        D = (imreconstruct(D+1,mask) - D)*FD.cellsize;
                    case 'carve'

                end
                if verbose
                    disp([datestr(clock) ' -- Weights for graydist calculated'])
                end

                % Here we calculate the auxiliary topography. That is, the
                % cost surface seeded at socalled PreSillPixels, i.e. the
                % pixel immediately upstream to sill pixels.
                D(I) = inf;
                D = graydist(double(D),double(PreSillPixel),'q') + 1;
                D(I) = -inf;

                if verbose
                    disp([datestr(clock) ' -- Auxiliary topography in flats calculated'])
                end
                
                
                if mexx
                    % mexed version: 
                    % identifies steepest downward neighbors in the DEM and 
                    % the distance grid obtained from graydist and performs
                    % a topological sort
                    clear I S demf                     
                    D(isinf(D)) = inf;
                    D(SILLS)    = 0;
                    clear SILLS
                    dem = double(dem);
                    % call to steepest neighbor. 
                    SE = steepestneighbor_mex(dem,D,FD.cellsize);
                    if verbose
                       disp([datestr(clock) ' -- Steepest neighbor identified'])
                    end
                    clear D
                    I  = isnan(dem);
                    clear dem
                    IX = find(SE==0 & ~I);
                    MV = nnz(I);
                    [FD.ix,FD.ixc] = tsort_mex(SE,uint32(IX),uint32(MV));
                    I  = FD.ix==0;
                    FD.ix(I) = [];
                    FD.ixc(I) = [];
                    
                else

                    % Find steepest neighbor
                    % sort pixels in D
                    % adapted from sortrows.m
                    clear PreSillPixel
                    [~,IXSortedFlats] = sort(D(:),'descend');
                    clear D

                    if verbose
                        disp([datestr(clock) ' -- Pixels sorted (1)'])
                    end

                    ndx = (uint32(1):uint32(nrc))';
                    ndx = ndx(IXSortedFlats);
                    clear IXSortedFlats

                    [~,FD.ix] = sort(dem(ndx),'descend');

                    if verbose
                        disp([datestr(clock) ' -- Pixels sorted (2)'])
                    end

                    FD.ix = uint32(FD.ix);
                    FD.ix = ndx(FD.ix);
                    clear ndx;



                    % a fast solution that has quite much memory overhead...
                    pp = zeros(FD.size,'uint32');
                    IX = (uint32(1):uint32(numel(dem)))';
                    pp(FD.ix) = IX;

                    % cardinal neighbors
                    IXC1 = imdilate(pp,[0 1 0; 1 1 1; 0 1 0]>0);
                    xxx1 = IXC1;
                    IX   = IXC1(FD.ix);
                    IXC1 = FD.ix(IX);
                    G1   = (dem(FD.ix)-dem(IXC1))/(FD.cellsize);
                    G1(FD.ix == IXC1) = -inf;

                    % diagonal neighbors
                    IXC2 = imdilate(pp,[1 0 1; 0 1 0; 1 0 1]>0);
                    xxx2 = IXC2;
                    IX   = IXC2(FD.ix);
                    IXC2 = FD.ix(IX);
                    G2   = (dem(FD.ix)-dem(IXC2))/(norm([FD.cellsize,FD.cellsize]));

                    % choose the steeper one

                    I    = G1<=G2 & xxx2(FD.ix)>xxx1(FD.ix);
                    FD.ixc = IXC1;
                    FD.ixc(I) = IXC2(I);

                    I = FD.ixc == FD.ix;
                    FD.ix(I) = [];
                    FD.ixc(I) = [];

                    % remove nans
                    I = isnan(dem);
                    FD.ixc(I(FD.ix)) = [];
                    FD.ix(I(FD.ix)) = [];
                end

                if verbose
                    disp([datestr(clock) ' -- Ordered topology established'])
                end
                
                        
            end
            
        else
            %% Multiple flow direction
            FD.cellsize = DEM.cellsize;
            FD.refmat   = DEM.refmat;
            FD.georef   = DEM.georef;
            FD.size     = DEM.size;
            
            switch lower(varargin{1});
                case 'multi'           
                    M = flowdir(DEM,'type','multi');                  
                    FD.type = 'multi';
                case 'dinf'
                    R = dem_flow(DEM.Z);
                    M = flow_matrix(DEM.Z,R,DEM.cellsize,DEM.cellsize);
                    M = -M';
                    
                    FD.type = 'Dinf';
                otherwise
                    error('unknown flow direction type')
            end
            
            try
                p = toposort(digraph(M));
                p = uint32(p);
            catch
                [p,~] = dmperm(speye(size(M))-M);
                p = uint32(p);
            end
            
            [FD.ix,FD.ixc,FD.fraction] = find(M(p,p));
            
            
            p = p(:);
            FD.ixc = p(FD.ixc);
            FD.ix  = p(FD.ix);
            
        end
    end


    function FD = saveobj(FD)
        FD.fastindexing = false;
    end
    
    function FD = set.fastindexing(FD,val)
        % fastindexing enables quick traversal along flow directions
        % starting from a single pixel. It is used when deriving flow 
        % paths, e.g., flowpathextract.
        
        validateattributes(val,{'numeric','logical'},{'scalar'})
        
        if val
            if ~strcmp(FD.type,'single');
                error('TopoToolbox:fastindexing','Fast indexing is only possible for single flow directions');
            end
            FD.fastindexing = true;
            FD.ixcix  = zeros(FD.size,'uint32');
            FD.ixcix(FD.ix) = uint32(1):uint32(numel(FD.ix));
        else
            FD.fastindexing = false;
            FD.ixcix  = [];
        end
    end
    
    function FD = multi_normalize(FD)
        
        if isempty(FD.fraction)
            return
        end
        s = accumarray(FD.ix,FD.fraction,[prod(FD.size) 1],@sum);
        FD.fraction = FD.fraction./s(FD.ix);
    end
    
    function FD = multi_weights(FD,varargin)
        
        p = inputParser;
        addParameter(p,'DEM',[],@(x) isa(x,'GRIDobj'));
        addParameter(p,'type','random')
        addParameter(p,'beta',1)
        parse(p,varargin{:});
        
        switch p.Results.type
            case 'random'
                FD.fraction = rand(size(FD.fraction));
        end
        FD = multi_normalize(FD);
        
        
        
    end
    
            

end
end


