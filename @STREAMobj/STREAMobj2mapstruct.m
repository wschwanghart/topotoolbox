function [GS,x,y] = STREAMobj2mapstruct(S,varargin)

%STREAMOBJ2MAPSTRUCT convert instance of STREAMobj to mapstruct
%
% Syntax
%
%     MS = STREAMobj2mapstruct(S)
%     MS = STREAMobj2mapstruct(S,type)
%     MS = STREAMobj2mapstruct(S,'seglength',length,...
%                        'attributes',{'fieldname1' var1 aggfunction1 ...
%                                      'fieldname2' var2 aggfunction2})
%     [MS,x,y] = ...
%     
%
% Description
%
%     A mapstruct is a structure array that contains vector geographic
%     data. STREAMobj2mapstruct converts an instance of STREAMobj to a
%     mapstruct MS. MS can be exported to a shapefile using the function
%     shapewrite available with the Mapping Toolbox. It can be plotted with
%     the function mapshow.
%
%     When called with following syntax
%     MS = STREAMobj2mapstruct(S)
%     then
%     MS contains individual line features for each river reach with unique
%     streamorder and contains following fields:
%
%        Geometry: 'Line'
%               X: vector with vertices of x-coordinates
%               Y: vector with vertices of y-coordinates
%     streamorder: stream order (Strahler)
%              IX: stream id
%        tribtoIX: tributary to stream with id IX
%
%     When called with following syntax
%     MS = STREAMobj2mapstruct(S,'seglength',length,...
%                        'attributes',{'fieldname1' var1 aggfunction1 ...
%                                      'fieldname2' var2 aggfunction2})
%
%     the streamnetwork is subdivided in reaches with approximate length in
%     map units defined by the parameter value pair 'seglength'-length. In
%     addition, a number of additional variables (instance of GRIDobjs or
%     node attributes (e.g. as returned by STREAMobj/streamorder,
%     STREAMobj/gradient, etc.) can be defined to be written as attributes
%     into the structure array. Since the values associated with one line
%     segment must be scalars, an aggregation function must be defined that
%     takes a vector and returns a scalar (e.g. @mean, @max, @std).
%
%     Note that the segment length in the output structure array is
%     variable since the function tries to evenly distribute the segments
%     between confluences, confluences and outlets, and confluences and
%     channelheads.
%
%     The additional output arguments x and y are nan-separated vectors
%     that can be used to plot the line features. Note that the order of x
%     and y is different from the output of STREAMobj2XY. 
%
% Input arguments
%
%     S       STREAMobj
%     type    {'strahler'} or 'shreve'
%
% Parameter name/value pairs
%
%     'seglength'   approximate segment length (scalar, map units)
%     'attributes'  1 x n*3 cell array where n is the number of attributes,
%                   and elements refer repeatively to
%                   {'fieldname' var aggfunction ...}
%                   where fieldname is the field name in the attribute
%                   table
%                   var is a GRIDobj associated with the STREAMobj or a
%                   node attribute list
%                   aggfunction is a anonymous function that takes a vector
%                   and returns a scalar
%                   *** This might not be easy to understand so check the 
%                       example below ***
%      'parallel'   {false} or true. Requires the Parallel Computing
%                   Toolbox. If set to true, the function will split the
%                   stream network into its connected components and
%                   process each individually in parallel. Mostly, there
%                   will not be a gain in speed as there is a significant
%                   computational overhead to split the network. Running in 
%                   may be more efficient if network is very large with 
%                   numerous connected components.   
%      'latlon'     {false} or true. If true, the function will attempt to
%                   transform coordinates to geographic coordinates. This
%                   works only if S contains a valid projection. 
%     
% Output arguments
%
%     MS      mapstruct
%     x       x coordinate vector
%     y       y coordinate vector
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     A = flowacc(FD);
%     S = STREAMobj(FD,A>1000);
%     DEM = imposemin(FD,DEM);
%     g = gradient(S,DEM);
%     MS = STREAMobj2mapstruct(S,'seglength',300,...
%                 'attributes',{'flowacc' A @mean 'gradient' g @mean});
%
%     % then export MS as shapefile using shapewrite
%     
% See also: shapewrite, STREAMobj2XY
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 27. September, 2021




if nargin == 1
    type = 'strahler';
elseif nargin == 2
    type = validatestring(varargin{1},{'strahler','shreve'},'STREAMobj/streamorder','type',2);
end

if nargin <= 2
    % get streamorder
    s = streamorder(S,type);
    % get outlets
    outlets = streampoi(S,'outlets','logical');
    
    % channel index counter
    chc = 1;
    CHN = zeros(size(outlets));
    
    % create empty mapstruct
    GS = struct('Geometry',{},...
        'X',{},...
        'Y',{},...
        'streamorder',{},...
        'IX',{},...
        'tribtoIX',{});
    
    
    for r = numel(S.ix):-1:1
        
        if CHN(S.ixc(r)) ==  0
            CHN(S.ixc(r)) = chc;
            GS(CHN(S.ixc(r))).Geometry = 'Line';
            GS(CHN(S.ixc(r))).X = S.x(S.ixc(r));
            GS(CHN(S.ixc(r))).Y = S.y(S.ixc(r));
            GS(CHN(S.ixc(r))).streamorder = s(S.ixc(r));
            GS(CHN(S.ixc(r))).IX = chc;
            
            chc = chc+1;
        end
        
        if s(S.ixc(r)) == s(S.ix(r))
            CHN(S.ix(r)) = CHN(S.ixc(r));
            
            GS(CHN(S.ix(r))).X(end+1) = S.x(S.ix(r));
            GS(CHN(S.ix(r))).Y(end+1) = S.y(S.ix(r));
            
            CHN(S.ixc(r)) = CHN(S.ix(r));
            
        elseif s(S.ixc(r)) ~= s(S.ix(r))
            
            
            CHN(S.ix(r)) = chc;
            GS(CHN(S.ix(r))).Geometry = 'Line';
            GS(CHN(S.ix(r))).X = S.x(S.ixc(r));
            GS(CHN(S.ix(r))).Y = S.y(S.ixc(r));
            
            GS(CHN(S.ix(r))).X(end+1) = S.x(S.ix(r));
            GS(CHN(S.ix(r))).Y(end+1) = S.y(S.ix(r));
            
            GS(CHN(S.ix(r))).IX = CHN(S.ix(r));
            GS(CHN(S.ix(r))).streamorder = s(S.ix(r));
            GS(CHN(S.ix(r))).tribtoIX = CHN(S.ixc(r));
            
            chc = chc+1;
        end
    end
    
    for r = 1:numel(GS)
        GS(r).X = [GS(r).X(end:-1:1) nan];
        GS(r).Y = [GS(r).Y(end:-1:1) nan];
        if isempty(GS(r).tribtoIX)
            GS(r).tribtoIX = 0;
        end
        
    end
else
    
    % do the first series of input argument checking
    p = inputParser;
    p.FunctionName = 'STREAMobj2mapstruct';
    addRequired(p,'S',@(x) isa(x,'STREAMobj'));
    addParameter(p,'seglength',S.cellsize*10,@(x) isscalar(x) && x>=S.cellsize*3);
    addParameter(p,'attributes',{},@(x) iscell(x)); 
    addParameter(p,'checkattributes',true,@(x) isscalar(x))
    addParameter(p,'parallel',false)
    addParameter(p,'latlon',false)
    parse(p,S,varargin{:});
    
    % check the attributes
    attributes   = p.Results.attributes;
    nrattributes = numel(attributes)/3;
    runinpar     = p.Results.parallel;
    seglength    = p.Results.seglength;
    
    % make valid and unique field names if required
    if p.Results.checkattributes
        attributes(1:3:end) = matlab.lang.makeValidName(attributes(1:3:end));
        attributes(1:3:end) = matlab.lang.makeUniqueStrings(attributes(1:3:end));
    end
    
    % convert strings to functions, if necessary
    for r = 1:nrattributes
        ix = ((r-1)*3)+3;
        if ischar(attributes{ix})
            attributes{ix} = str2func(attributes{ix});
        end
    end
    
    % convert attributes to node attribute lists
    for r = 1:nrattributes
        ix = ((r-1)*3)+2;
        if isa(attributes{ix},'GRIDobj')
            attributes{ix} = getnal(S,attributes{ix});
        elseif isnal(S,attributes{ix})
            % everything's fine here
        else
            error('TopoToolbox:incompatibleformat',...
                ['Incompatible format of attribute. Attribute is expected to\n'...
                 'be either a GRIDobj or a node attribute list (nal) of S']);
        end
    end
    
    %% Parallel computing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if runinpar
        % Place each connected component in an individual cell
        [CS,locS] = STREAMobj2cell(S);
        
        % Derive a cell array of attributes
        CATT = cell(numel(CS),nrattributes*3);
        for r = 1:numel(CS)
            for r2 = 1:nrattributes
                CATT{r,(r2-1)*3 + 1} = attributes{(r2-1)*3 + 1}; 
                CATT{r,(r2-1)*3 + 3} = attributes{(r2-1)*3 + 3};               
                CATT{r,(r2-1)*3 + 2} = attributes{(r2-1)*3 + 2}(locS{r});
            end
        end
        
        % Preallocate cell array
        GS = cell(numel(CS),1);
        
        parfor r = 1:numel(CS)            
            GS{r} = STREAMobj2mapstruct(CS{r},'seglength',seglength,...
                                                'attributes',CATT(r,:),...
                                                'parallel',false);
        end
        
        GS = vertcat(GS{:});
        if nargout>1
            x = [GS.X]';
            y = [GS.Y]';
        end
        return
        
    end
    %% Parallel computing ends here %%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            
    
    % write stream segments
    confluences = streampoi(S,'confl','logical');
    
    % preallocate cell array with attributes
    attr = cell(1,nrattributes);
    
    % get nan-punctuated vectors
    [x,y,confluences,attr{:}] = STREAMobj2XY(S,confluences,attributes{2:3:end});
    
    % in addition to existing nan-separators, separate vectors at confluences
    % the confluences which are followed by nan can be neglegted
    confluences(circshift(isnan(confluences),-1) & confluences==1) = 0;
    
    % doublicate nodes at remaining confluences (insert rows)
    ind = find(confluences == 1);
    x  = insertrows(x,x(ind),ind);
    y  = insertrows(y,y(ind),ind);
    confluences = insertrows(confluences,0,ind);
    for r = 1:numel(attr)
        attr{r} = insertrows(attr{r},attr{r}(ind),ind);
    end
    
    % insert nans
    ind = find(confluences == 1);
    x  = insertrows(x,nan,ind);
    y  = insertrows(y,nan,ind);
    for r = 1:numel(attr)
        attr{r} = insertrows(attr{r},nan,ind);
    end
    
    % get additional separators to segment lines further according to
    % segment length
    % calculate the distance within each individual segment
    sep = isnan(x);
    d  = zeros(size(x));
    for r = 2:numel(d)
        if sep(r-1)
            d(r) = 0;
        elseif sep(r)
            d(r) = nan;
        else
            d(r) = d(r-1)+sqrt((x(r)-x(r-1)).^2 + (y(r)-y(r-1)).^2);
        end
    end
    
    % try to get equal length segments being equally distributed over an
    % existing segment.
    % Use equal quantiles
    ixsep = find(sep);
    nrsep = numel(ixsep);
    ixsep = [0;ixsep];
    newsep   = zeros(size(d));
    for r = 1:nrsep
        dd = d(ixsep(r)+1 : ixsep(r+1)-1);
        ddmax = max(dd);
        nrnewsegs = max(round(ddmax/p.Results.seglength),1);
        
        if nrnewsegs > 1
            [~,segtemp] = histc(dd,linspace(0,ddmax + S.cellsize*.01,nrnewsegs));
            segtemp = [0;diff(segtemp)];
            newsep(ixsep(r)+1 : ixsep(r+1)-1) = segtemp;
        end
        
    end
    
    % now set new separators according to newsep
    ind = find(newsep > 0);
    x  = insertrows(x,x(ind),ind);
    y  = insertrows(y,y(ind),ind);
    newsep = insertrows(newsep,0,ind);
    for r = 1:numel(attr)
        attr{r} = insertrows(attr{r},attr{r}(ind),ind);
    end
    
    % insert nans
    ind = find(newsep > 0);
    x  = insertrows(x,nan,ind);
    y  = insertrows(y,nan,ind);
    for r = 1:numel(attr)
        attr{r} = insertrows(attr{r},nan,ind);
    end
    
    
    % Now, finally, create the mapstruct
    ix      = find(isnan(x));
    nrlines = numel(ix);
    % make attributes a row (cell) vector
    attributes = attributes(:)';
    cc  = [attributes(1:3:end);
           repmat({cell(nrlines,1)},1,nrattributes)];
    
    GS = struct('Geometry','Line',...
                'X',cell(nrlines,1),...
                'Y',cell(nrlines,1),...    
                cc{:});
    
    IXs = [1;ix(1:end-1)+1];
    IXe = ix-1;
    
    
    for r = 1:nrlines
        GS(r).X = x(IXs(r):IXe(r));
        GS(r).Y = y(IXs(r):IXe(r));
        
        for r2 = 1:nrattributes
            GS(r).(attributes{((r2-1)*3)+1}) = attributes{((r2-1)*3)+3}(attr{r2}(IXs(r):IXe(r)-1));
        end
        
    end
    
end
    
% Convert to lat lon if requested
if exist('p', 'var')
if p.Results.latlon
    for r = 1:numel(GS)
        [GS(r).Y,GS(r).X] = minvtran(S.georef.mstruct,GS(r).X,GS(r).Y);
    end
end
end

if nargout>1
    x = [GS.X]';
    y = [GS.Y]';
end




end



%% INSERTROWS by JOS
% http://www.mathworks.com/matlabcentral/fileexchange/9984-insertrows-v2-0-may-2008
function [C,RA,RB] = insertrows(A,B,ind)
% INSERTROWS - Insert rows into a matrix at specific locations
%   C = INSERTROWS(A,B,IND) inserts the rows of matrix B into the matrix A at
%   the positions IND. Row k of matrix B will be inserted after position IND(k)
%   in the matrix A. If A is a N-by-X matrix and B is a M-by-X matrix, C will
%   be a (N+M)-by-X matrix. IND can contain non-integers.
%
%   If B is a 1-by-N matrix, B will be inserted for each insertion position
%   specified by IND. If IND is a single value, the whole matrix B will be
%   inserted at that position. If B is a single value, B is expanded to a row
%   vector. In all other cases, the number of elements in IND should be equal to
%   the number of rows in B, and the number of columns, planes etc should be the
%   same for both matrices A and B. 
%
%   Values of IND smaller than one will cause the corresponding rows to be
%   inserted in front of A. C = INSERTROWS(A,B) will simply append B to A.
%
%   If any of the inputs are empty, C will return A. If A is sparse, C will
%   be sparse as well. 
%
%   [C, RA, RB] = INSERTROWS(...) will return the row indices RA and RB for
%   which C corresponds to the rows of either A and B.
%
%   Examples:
%     % the size of A,B, and IND all match
%        C = insertrows(rand(5,2),zeros(2,2),[1.5 3]) 
%     % the row vector B is inserted twice
%        C = insertrows(ones(4,3),1:3,[1 Inf]) 
%     % matrix B is expanded to a row vector and inserted twice (as in 2)
%        C = insertrows(ones(5,3),999,[2 4])
%     % the whole matrix B is inserted once
%        C = insertrows(ones(5,3),zeros(2,3),2)
%     % additional output arguments
%        [c,ra,rb] = insertrows([1:4].',99,[0 3]) 
%        c.'     % -> [99 1 2 3 99 4] 
%        c(ra).' % -> [1 2 3 4] 
%        c(rb).' % -> [99 99] 
%
%   Using permute (or transpose) INSERTROWS can easily function to insert
%   columns, planes, etc:
%
%     % inserting columns, by using the transpose operator:
%        A = zeros(2,3) ; B = ones(2,4) ;
%        c = insertrows(A.', B.',[0 2 3 3]).'  % insert columns
%     % inserting other dimensions, by using permute:
%        A = ones(4,3,3) ; B = zeros(4,3,1) ; 
%        % set the dimension on which to operate in front
%        C = insertrows(permute(A,[3 1 2]), permute(B,[3 1 2]),1) ;
%        C = ipermute(C,[3 1 2]) 
%
%  See also HORZCAT, RESHAPE, CAT

% for Matlab R13
% version 2.0 (may 2008)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History:
% 1.0, feb 2006 - created
% 2.0, may 2008 - incorporated some improvements after being selected as
% "Pick of the Week" by Jiro Doke, and reviews by Tim Davis & Brett:
%  - horizontal concatenation when two arguments are provided
%  - added example of how to insert columns
%  - mention behavior of sparse inputs
%  - changed "if nargout" to "if nargout>1" so that additional outputs are
%    only calculated when requested for

narginchk(2,3);

if nargin==2,
    % just horizontal concatenation, suggested by Tim Davis
    ind = size(A,1) ;
end

% shortcut when any of the inputs are empty
if isempty(B) || isempty(ind),    
    C = A ;     
    if nargout > 1,
        RA = 1:size(A,1) ;
        RB = [] ;
    end
    return
end

sa = size(A) ;

% match the sizes of A, B
if numel(B)==1,
    % B has a single argument, expand to match A
    sb = [1 sa(2:end)] ;
    B = repmat(B,sb) ;
else
    % otherwise check for dimension errors
    if ndims(A) ~= ndims(B),
        error('insertrows:DimensionMismatch', ...
            'Both input matrices should have the same number of dimensions.') ;
    end
    sb = size(B) ;
    if ~all(sa(2:end) == sb(2:end)),
        error('insertrows:DimensionMismatch', ...
            'Both input matrices should have the same number of columns (and planes, etc).') ;
    end
end

ind = ind(:) ; % make as row vector
ni = length(ind) ;

% match the sizes of B and IND
if ni ~= sb(1),
    if ni==1 && sb(1) > 1,
        % expand IND
        ind = repmat(ind,sb(1),1) ;
    elseif (ni > 1) && (sb(1)==1),
        % expand B
        B = repmat(B,ni,1) ;
    else
        error('insertrows:InputMismatch',...
            'The number of rows to insert should equal the number of insertion positions.') ;
    end
end

sb = size(B) ;

% the actual work
% 1. concatenate matrices
C = [A ; B] ;
% 2. sort the respective indices, the first output of sort is ignored (by
% giving it the same name as the second output, one avoids an extra 
% large variable in memory)
[~,abi] = sort([[1:sa(1)].' ; ind(:)]) ;
% 3. reshuffle the large matrix
C = C(abi,:) ;
% 4. reshape as A for nd matrices (nd>2)
if ismatrix(A),
    sc = sa ;
    sc(1) = sc(1)+sb(1) ;
    C = reshape(C,sc) ;
end

if nargout > 1,
    % additional outputs required
    R = [zeros(sa(1),1) ; ones(sb(1),1)] ;
    R = R(abi) ;
    RA = find(R==0) ;
    RB = find(R==1) ;
end
end










