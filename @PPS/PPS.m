classdef PPS

%PPS Point patterns on stream networks
%
% Syntax
%
%     P = PPS(S,'PP',IX,...)
%     P = PPS(S,'PP',xy,...)
%     P = PPS(S,'PP',MS,...)
%     P = PPS(S,'PP',nal,...)
%     P = PPS(S,'PP',geotable,...)
%     P = PPS(S,'runif',n,...)
%     P = PPS(S,'rpois',lambda,...)
%     P = PPS(S,'intersect',PS,...)
%     P = PPS(...,'z',DEM)
%
% Description
%
%     An instance of class PPS stores a stream network (STREAMobj) and an
%     associated point pattern. PPS is also the constructor function.
%     There are following syntaxes to create an instance of PPS:
%
%     P = PPS(S,'PP',IX,...) creates a point pattern from points with the 
%     linear index IX into a GRIDobj (see TopoToolbox).
%
%     P = PPS(S,'PP',xy,...) uses the coordinates in the px2 matrix with p
%     points and x-coordinates in the first column and y-coordinates in the
%     second column. Coordinates that are not located on the stream network
%     are snapped using STREAMobj/snap2stream(S,...). 
%
%     P = PPS(S,'PP',MS,...) takes the mapping structure (as returned by
%     shaperead). The mapping structure must be a point vector shape with
%     the same coordinate system as S. Coordinates will be snapped to
%     streams.
%
%     P = PPS(S,'PP',nal,...) extracts points from the logical node-
%     attribute list nal. The vector nal must have as many elements as
%     there are nodes in stream network.
%
%     P = PPS(S,'runif',n,...) creates a spatially uniformly distributed 
%     point pattern of n points on the stream network S.
%
%     P = PPS(S,'rpois',lambda,...) creates a point pattern on S which is
%     completely spatial random (CSR) (Poisson distribution) with intensity
%     lambda. Lambda is the average density of points (points per unit
%     length). 
%
%     P = PPS(S,'intersect',PS,...) creates a point pattern on S where
%     points are derived from the spatial intersections of the stream
%     network and lines or polygons in PS. PS can be a two-column matrix
%     with x and y coordinates (eventually separated by nans), a mapping
%     structure array, or a mapshape or geoshape object.
%
%     P = PPS(...,'z',DEM) adds elevations (either a GRIDobj or
%     node-attribute list compatible with S) to each vertex in the stream
%     network S. Elevations are necessary for functions such as plotdz(P)
%     to run properly.
%
% Input arguments
%
%     S       STREAMobj
%     IX      linear index into digital elevation model from which S was
%             derived.
%     xy      nx2 matrix of coordinates. Coordinates not on the stream
%             network will be snapped to the stream network (see function
%             STREAMobj/snap2stream).
%     nal     logical node-attribute list
%     n       number of points if binomial point process
%     lambda  point pattern intensity (i.e., points per m along the
%             network) if Poisson point process
%
%     Parameter name/value pairs
%
%     'z'            GRIDobj or node-attribute list of elevations
%     'alongstream'  only applicable if called together with 'PP',xy. 
%                    By default, []. If FLOWobj is supplied, then stream 
%                    snapping moves the points along the flow network (see
%                    STREAMobj/snap2stream).
%
% Output arguments
%
%     P       instance of class PPS
%
% Example 1: Create Poisson-distributed point pattern on stream network  
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     S = removeshortstreams(S,100);
%     S = clean(S);
%     P = PPS(S,'rpois',0.001,'z',DEM);
%     plot(P)
%
% Example 2: 
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S,1);
%
%     % Find knickpoints
%     [~,kp] = knickpointfinder(S,DEM,'tol',30,...
%                               'split',false,...
%                               'verbose',false,...
%                               'plot',false);
%     P = PPS(S,'PP',kp.IXgrid,'z',DEM);
%     plot(P)
% 
% Reference: Schwanghart, W., Molkenthin, C., & Scherler, D. (2020). A 
% systematic approach and software for the analysis of point patterns on 
% river networks. Earth Surface Processes and Landforms, 46 (9), 1847-1862. 
% [DOI: 10.1002/esp.5127]
%
% See also: GRIDobj, FLOWobj, STREAMobj, STREAMobj/snap2stream, SWATHobj
%           DIVIDEobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 21. December, 2021
    
properties
%S instance of STREAMobj
%    The S property contains a stream network stored as STREAMobj.
%
%    See also STREAMobj
    S   = STREAMobj
    
%PP column vector of doubles
%    PP contains the linear indices into the node-attribute list of S.
%
%    See also STREAMobj
    PP  = []
end

properties (Dependent = true) 
%z node-attribute list with elevation values
%    z contains a node-attribute list with elevation values
%
%    See also STREAMobj
    z
%covariates  table with covariates (still unsupported)
    covariates
%marks   table with marks (still unsupported)    
    marks
%ppxy   dynamic property with x and y coordinates of points
    ppxy
%ppz    point elevations
    ppz
%ppcovariates   covariate values at point locations (still unsupported)
    ppcovariates
end

properties (Access=private)
    Private_z = []
    Private_covariates = table
    Private_marks = table
end


methods
    function [P,II] = PPS(S,varargin)
        
        narginchk(1,inf);
               
        % Check input arguments
        p = inputParser;
        p.FunctionName = 'PPS';
        addRequired(p,'S',@(x) isa(x,'STREAMobj'));
        
        % Add elevation
        addParameter(p,'z',[],@(x) isa(x,'GRIDobj') || isnal(S,x));
        
        % Either, points are provided (this may include marks, but this is not yet supported)
        addParameter(p,'PP',[]);
        addParameter(p,'marks',table, @(x) istable(x));
        
        % Or the users sets a number of options to simulate points
        % 1. Number of uniform points
        addParameter(p,'runif',[],@(x) isscalar(x)); 
        % 2. Intensity of poisson process
        addParameter(p,'rpois',[],@(x) isnumeric(x));
        
        % Or find intersections with polylines or polygons
        addParameter(p,'intersect',[]);
        
        % Finally, covariates may be included
        addParameter(p,'covariates',table,@(x) istable(x));
        addParameter(p,'interactive',false);
        
        % Some additional parameters
        addParameter(p,'warning',true);
        addParameter(p,'alongflow',[],@(x) isa(x,'FLOWobj'));

        % Parse
        parse(p,S,varargin{:});
        
        results = p.Results;
        P.S = results.S;
        
        %% Elevations
        if isa(results.z,'GRIDobj')
            P.z = getnal(S,results.z);
        else
            P.z = results.z;
        end

        %% Enable reading of geotable
        if ~verLessThan('map','5.2')
            if isgeotable(results.PP)
                points = results.PP.Shape;
                x = [points.X];
                y = [points.Y];

                results.PP = [x(:) y(:)];
            end
        end
        
        %% Read or create points
        if ~isempty(results.PP) 
            %% Points (PP)
            % Users can submit different kind of points
            if iscolumn(results.PP) && ~islogical(results.PP) && ~isstruct(results.PP)
                % linear index
                IXgrid     = results.PP;
                [II,P.PP]   = ismember(IXgrid,P.S.IXgrid);
                if ~all(II) && p.Results.warning
                    warning([num2str(numel(II)-nnz(II)) ' of ' num2str(numel(II)) ...
                        ' points are not on the stream network'])
                    P.PP = P.PP(II);
                end

            elseif isnal(S,results.PP) && islogical(results.PP)
                % logical node attribute list
                P.PP = find(results.PP);
                
            elseif (size(results.PP,2) == 2 && ~islogical(results.PP)) || isstruct(results.PP)
                % structure array or array with coordinates.
                if isstruct(results.PP)
                    x = [results.PP.X]';
                    y = [results.PP.Y]';
                else
                    x = results.PP(:,1);
                    y = results.PP(:,2);
                end
                
                % coordinates
                [~,~,IXgrid] = snap2stream(S,x,y,'alongflow',results.alongflow);
                [~,P.PP] = ismember(IXgrid,P.S.IXgrid);
            end
            
        %% Generate points            
        elseif ~isempty(results.runif)
            %% Uniform random
            IX = randlocs(P.S,results.runif,false);
            [~,P.PP] = ismember(IX,P.S.IXgrid);
        elseif ~isempty(results.rpois)
            %% Poisson distributed
            IX = randpoi(P.S,results.rpois);
            [~,P.PP] = ismember(IX,P.S.IXgrid);
        elseif ~isempty(results.intersect)
            %% From intersection with other data
            PS = results.intersect;
            if isnumeric(PS)
                X1 = PS(:,1);
                Y1 = PS(:,2);
            elseif isa(PS,'mapshape')
                X1 = PS.X;
                Y1 = PS.Y;
            elseif isa(PS,'struct')
                X1 = [PS.X];
                Y1 = [PS.Y];
            elseif isa(PS,'geoshape')
                lat = PS.Latitude;
                lon = PS.Longitude;
                [X1,Y1] = mfwdtran(P.S.georef.mstruct,lat,lon);
            end
            [X,Y] = STREAMobj2XY(P.S);
            [Xi,Yi] = polyxpoly(X1,Y1,X,Y);
            
            [~,~,IX] = snap2stream(P.S,Xi,Yi,'maxdist',P.S.cellsize,...
                'plot',false);
            [~,P.PP] = ismember(IX,S.IXgrid); 
        end
               
        %% Additional input argument checking for covariates and marks
        notallowed_varn = {'x','y','z'};
        % Marks (only applicable if PP is set, table must have the same 
        % height as PP)
        
        if ~isempty(results.marks)
            % make sure that there are as many observations as there are
            % points
            t = results.marks;
            
            if exist('II','var')
                t = t(II,:);
            end
            
            if height(t) ~= numel(P.PP)
                error('The number of points and marks must be the same.')
            end
            % some variable names should not be set in the table            
            for r = 1:numel(notallowed_varn)
                if any(strcmp(t.Properties.VariableNames,notallowed_varn{r}))
                    error(['Names of mark must not be ' notallowed_varn{r} '.']) 
                end
            end            
            P.Private_marks = t;
        end
               
        % Covariates (must be node attribute lists)
        if ~isempty(results.covariates)
            % make sure that there are as many observations as there are
            % points
            t = results.covariates;
            if height(t) ~= numel(P.S.IXgrid)
                error('Covariate table must have as many rows as there are nodes in the stream network.')
            end
            % some variable names should not be set in the table            
            for r = 1:numel(notallowed_varn)
                if any(strcmp(t.Properties.VariableNames,notallowed_varn{r}))
                    error(['Names of covariates must not be ' notallowed_varn{r} '.']) 
                end
            end            
            P.Private_covariates = t;
        end
        
    end
    
%     function P = subsref(P,S)
%         if numel(S) == 1
%         switch S(1).type
%             case '()'
%                 if numel(S.subs) == 1
%                     ix = S.subs{1};
%                     ix = ix(:);
%                     P.PP = P.PP(ix);
%                 else
%                     error('Subscript indexing not allowed.')
%                 end
%         end
%             case '.'
%                 P = P.(S.subs);
%             otherwise
%                 error('Only linear indexing with () allowed')
%         end
%     end
    
    %% methods to access and set properties
    function P = set.z(P,DEM)
        
        ST = P.S;
        
        if isa(DEM,'GRIDobj')
            P.Private_z = double(getnal(ST,DEM));
        else
            if isempty(DEM)
                P.Private_z = [];
                return
            end
            
            if ~isnal(ST,DEM)
                error('DEM must be either a DEM or a node-attribute list of elevation values')
            end
            P.Private_z = double(DEM);
        end
        
    end
    
    function val = get.z(P)       
        val = P.Private_z;
    end
    
    function value = get.ppxy(P)
        value = [P.S.x(P.PP) P.S.y(P.PP)];
    end
        
    function value = get.ppz(P)
        if isempty(P.z)
            value = [];
        else
            value = P.z(P.PP);
        end
    end    
    
    function value = get.ppcovariates(P)
        if isempty(P.covariates)
            value = table;
        else
            value = P.covariates(P.PP,:);
        end
    end
        
    function value = get.covariates(P)
        value = P.Private_covariates;
    end
           
    
end
end


