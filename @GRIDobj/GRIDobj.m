classdef GRIDobj
    
%GRIDobj Create instance of a GRIDobj
%
% Syntax
%
%     DEM = GRIDobj(X,Y,dem)
%     DEM = GRIDobj('ESRIasciiGrid.txt')
%     DEM = GRIDobj('GeoTiff.tif')
%     DEM = GRIDobj();
%     DEM = GRIDobj([]);
%     DEM = GRIDobj(FLOWobj or GRIDobj or STREAMobj,class)
%
%
% Description
%
%     GRIDobj creates an instance of the grid class, which contains a
%     numerical or logical matrix and information on georeferencing. When a
%     GRIDobj is created from a file, the number format of the data in
%     GRIDobj is either single or double. Unsigned and signed integers are
%     converted to single. For unsigned integers, missing values are
%     assumed to be denoted as intmax(class(input)). For signed integers,
%     missing values are assumed to be intmin(class(input)). Please check,
%     that missing values in your data have been identified correctly
%     before further analysis.
%
%     Note that while throughout this help text GRIDobj is associated with
%     gridded digital elevation models, instances of GRIDobj can contain
%     other gridded, single band, datasets such as flow accumulation grids, 
%     gradient grids etc.
%
%     DEM = GRIDobj(X,Y,dem) creates a DEM object from the coordnate
%     matrices or vectors X and Y and the matrix dem. The elements of dem
%     refer to the elevation of each pixel. 
%
%     DEM = GRIDobj('ESRIasciiGrid.txt') creates a DEM object from an ESRI 
%     Ascii grid exported from other GI systems. 
%
%     DEM = GRIDobj('GeoTiff.tif') creates a DEM object from a Geotiff.
%
%     DEM = GRIDobj() opens a dialog box to read either an ESRI Ascii Grid
%     or a Geotiff.
%
%     DEM = GRIDobj([]) creates an empty instance of GRIDobj
%
%     DEM = GRIDobj(FLOWobj or GRIDobj or STREAMobj,class) creates an
%     instance of GRIDobj with all common properties (e.g., spatial
%     referencing) inherited from another instance of a FLOWobj, GRIDobj 
%     or STREAMobj class. DEM.Z is set to all zeros where class can be
%     integer classes or double or single. By default, class is double.
%
% Example
%
%     % Load DEM
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     % Display DEM
%     imageschs(DEM)
%
% See also: FLOWobj, STREAMobj, GRIDobj/info
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017

    
    properties
        %Public properties
        
%Z matrix with elevation values
%    The Z property contains the elevation values in a 2D matrix.
%
%    See also GRIDobj
        Z         
        
%CELLSIZE cellsize of the grid (scalar)
%    The cellsize property specifies the spacing of the grid in x and y
%    directions. Note that TopoToolbox requires grids to have square cells,
%    e.g., dx and dy are the same.
%
%    See also GRIDobj       
        cellsize  
        
%REFMAT 3-by-2 affine transformation matrix
%    The refmat property specifies a 3-by-2 affine transformation matrix as
%    used by the mapping toolbox. 
%
%    See also GRIDobj, makerefmat
        refmat  
        
%SIZE size of the grid (two element vector)
%    The cellsize property is a two element vector that contains the number
%    of rows and columns of the Z matrix.
%
%    See also GRIDobj, size
        size 
        
%NAME optional name (string)
%    The name property allows to specify a name of the grid. By default and
%    if the constructor is called with a filename, the name property is set
%    to the name of the file.
%
%    See also GRIDobj
        name     
        
%ZUNIT unit of grid values (string)
%    The zunit is optional and is used to store the physical unit (e.g. m)
%    of an instance of GRIDobj. This property is currently not fully
%    supported and TopoToolbox functions usually assume that the unit is in
%    meters and equals the xyunit property.
%
%    See also GRIDobj
        zunit     
        
%XYUNIT unit of the coordinates (string)
%    The xyunit is optional and is used to store the physical unit (e.g. m)
%    of the coordinates. This property is currently not fully
%    supported and TopoToolbox functions usually assume that the unit is in
%    meters and equals the zunit property.
%
%    See also GRIDobj            
        xyunit    
        
%GEOREF additional information on spatial referencing (structure array)
%    The georef property stores an instance of map.rasterref.MapCellsReference, 
%    a structure array GeoKeyDirectoryTag, and a mapping structure
%    (mstruct).
%
%    See also GRIDobj, geotiffinfo        
        georef    
        
    end
    
    methods
        function DEM = GRIDobj(varargin)
            % GRIDobj constructor
          
            if nargin == 3                
                %% GRIDobj is created from three matrices
                % GRIDobj(X,Y,dem)
                X = varargin{1};
                Y = varargin{2};
                
                if min(size(X)) > 1
                    X = X(1,:);
                end
                if min(size(Y)) > 1
                    Y = Y(:,1);
                end
                
                DEM.Z   = varargin{3};
                DEM.size = size(DEM.Z);
                
                if numel(X) ~= DEM.size(2) || numel(Y) ~= DEM.size(1)
                    error('TopoToolbox:GRIDobj',...
                          ['Coordinate matrices/vectors don''t fit the size of the \n'...
                           'the grid']);
                end
                
                if (Y(2)-Y(1)) > 0
                    % the y coordinate vector must be monotonically
                    % decreasing, so that the left upper edge of the DEM is
                    % north-west (on the northern hemisphere).
                    DEM.Z = flipud(DEM.Z);
                    Y = Y(end:-1:1);
                end
                
                dy = Y(2)-Y(1);
                dx = X(2)-X(1);
                
                if abs(abs(dx)-abs(dy))>1e-9
                    error('TopoToolbox:GRIDobj',...
                          'The resolution in x- and y-direction must be the same');
                end
                    
                
                DEM.refmat = double([0 dy;...
                              dx 0;...
                              X(1)-dx Y(1)-dy]);
                
                DEM.cellsize = dx;
                DEM.georef = [];
                DEM.name   = [];
            
                
            elseif nargin <= 2
                
                
                if nargin == 0
                    %% No input arguments. File dialog box will open and ask
                    % for a txt or tiff file as input
                    FilterSpec  = {'*.txt;*.asc;*.tif;*.tiff','supported file types (*.txt,*.asc,*.tif,*.tiff)';...
                                   '*.txt',   'ESRI ASCII grid (*.txt)';...
                                   '*.asc',   'ESRI ASCII grid (*.asc)';...
                                   '*.tif',   'GeoTiff (*.tif)';...
                                   '*.tiff',  'GeoTiff (*.tiff)';...
                                   '*.*',     'all files (*.*)'};
                        
                    DialogTitle = 'Select ESRI ASCII grid or GeoTiff';
                    [FileName,PathName] = uigetfile(FilterSpec,DialogTitle);
                
                    if FileName == 0
                        error('TopoToolbox:incorrectinput',...
                                'no file was selected')
                    end
                
                    filename = fullfile(PathName, FileName);
                    
                elseif nargin > 0
                    % One input argument
                    if isempty(varargin{1})
                        % if empty array than return empty GRIDobj
                        return
                    elseif isa(varargin{1},'GRIDobj') || ...
                           isa(varargin{1},'FLOWobj') || ...
                           isa(varargin{1},'STREAMobj') 
                        % empty GRIDobj
                        DEM = GRIDobj([]);
                        % find common properties of F and G and from F to G
                        pg = properties(DEM);
                        pf = properties(varargin{1});
                        p  = intersect(pg,pf);
                        for r = 1:numel(p)
                            DEM.(p{r}) = varargin{1}.(p{r});
                        end
                        if nargin == 1
                            cl = 'double';
                        else
                            cl = varargin{2};
                        end
                        
                        if strcmp(cl,'logical')
                            DEM.Z = false(DEM.size);
                        else
                            DEM.Z = zeros(DEM.size,cl);
                        end
                        DEM.name = '';
                            
                        return
                    end
                
                    % GRIDobj is created from a file
                    filename = varargin{1};
                end
                
                % check if file exists
                if exist(filename,'file')~=2
                    error('File doesn''t exist')
                end
                    
                
                % separate filename into path, name and extension
                [pathstr,DEM.name,ext]  = fileparts(filename);
                

                if any(strcmpi(ext,{'.tif', '.tiff'}))
                    % it is a GeoTiff
                    try 
                        % try to read using geotiffread (requires mapping
                        % toolbox)
                        [DEM.Z, DEM.refmat, ~] = geotiffread(filename);
                        gtiffinfo              = geotiffinfo(filename);
                        DEM.georef.SpatialRef  = gtiffinfo.SpatialRef; 
                        DEM.georef.GeoKeyDirectoryTag = gtiffinfo.GeoTIFFTags.GeoKeyDirectoryTag;
                        georef_enabled = true;
                        
                    catch ME
                        
                        % mapping toolbox is not available. Will try to
                        % read the tif file together with the tfw file
                        georef_enabled = false;
                        
                        % the tfw file has the same filename but a tfw
                        % extension
                        tfwfile = fullfile(pathstr,[DEM.name '.tfw']);
                        
                        % check whether file exists. If it exists then read
                        % it using worldfileread or own function
                        tfwfile_exists = exist(tfwfile,'file');
                        if tfwfile_exists
                            try 
                                % prefer builtin worldfileread, if
                                % available
                                DEM.refmat = worldfileread(tfwfile);
                            catch ME
                                W = dlmread(tfwfile);
                                
                                DEM.refmat(2,1) = W(1,1);
                                DEM.refmat(1,2) = W(4,1);
                                DEM.refmat(3,1) = W(5,1)-W(1);
                                DEM.refmat(3,2) = W(6,1)-W(4);
                            end
                            
                            DEM.Z = imread(filename);
                        else
                            if ~tfwfile_exists
                                error('TopoToolbox:GRIDobj:read',...
                                    'GRIDobj cannot read the TIF-file because it does not have a tfw-file.');
                            else
                                throw(ME)
                            end
                            
                        end
                    end
                    
                    
                    % Unless any error occurred, we now attempt to generate
                    % an mapping projection structure. This will not work
                    % if the DEM is in a geographic coordinate system or if
                    % the projection is not supported by mstruct.
                    if georef_enabled
                        try
                            DEM.georef.mstruct = geotiff2mstruct(gtiffinfo);
                        catch
                            DEM.georef.mstruct = [];
                            warning('TopoToolbox:GRIDobj:projection',...
                                ['GRIDobj cannot derive a map projection structure. This is either\n' ...
                                 'because the grid is in a geographic coordinate system or because\n' ...
                                 'geotiff2mstruct cannot identify the projected coordinate system used.\n' ...
                                 'TopoToolbox assumes that horizontal and vertical units of DEMs are \n'...
                                 'the same. It is recommended to use a projected coordinate system,\n' ...
                                 'preferably UTM WGS84. Use the function GRIDobj/reproject2utm\n' ...
                                 'to reproject your grid.'])
                        end
                    end
     
                    % Finally, check whether no_data tag is available. This tag is
                    % not accessible using geotiffinfo (nice hack by Simon
                    % Riedl)
                    tiffinfo = imfinfo(filename);
                    if isfield(tiffinfo,'GDAL_NODATA')
                        nodata_val = str2double(tiffinfo.GDAL_NODATA);
                    end
        
                else
                    [DEM.Z,R] = rasterread(filename);
                    DEM.refmat = R;
                    DEM.georef = [];
                end
                
                DEM.size = size(DEM.Z);
                DEM.cellsize = abs(DEM.refmat(2));
                
                % remove nans
                demclass = class(DEM.Z);
                nodata_val_exists = exist('nodata_val','var');
                
                switch demclass
                    case {'uint8','uint16','uint32'}
                        % unsigned integer
                        DEM.Z = single(DEM.Z);
                        
                        if nodata_val_exists
                            nodata_val = single(nodata_val);
                            DEM.Z(DEM.Z == nodata_val) = nan;
                        else                                                 
                            DEM.Z(DEM.Z==intmax(demclass)) = nan;
                        end
                        
                    case {'int8','int16','int32'}
                        % signed integer
                        DEM.Z = single(DEM.Z);
                        if nodata_val_exists
                            nodata_val = single(nodata_val);
                            DEM.Z(DEM.Z == nodata_val) = nan;
                        else                                                 
                            DEM.Z(DEM.Z==intmin(demclass)) = nan;
                        end
                        
                    case {'double','single'}
                        if nodata_val_exists
                            DEM.Z(DEM.Z == cast(nodata_val,class(DEM.Z))) = nan;
                        end
                    case 'logical'
                    otherwise
                        error('TopoToolbox:GRIDobj','unrecognized class')
                end
                
                
            end 
        end
    
    end
    
end
    



% Subfunction for ASCII GRID import
function [Z,refmat] = rasterread(file)

fid=fopen(file,'r');
% loop through header

header = struct('ncols',[],...
                'nrows',[],...
				'xllcorner',[],...
				'yllcorner',[],...
				'cellsize',[],...
				'nodata',[]);
				
names   = fieldnames(header);
nrnames = numel(names);

try
    fseek(fid,0,'bof');
    for r = 1:nrnames ;
        headertext = fgetl(fid);
        [headertext, headernum] = strtok(headertext,' ');
        I = cellfun(@(x,y) strcmpi(x(1:4),y(1:4)),names,repmat({headertext},nrnames,1));
        header.(names{I}) = str2double(headernum);
    end
catch ME1
    error('header can not be read')
end


% read raster data
Z = fscanf(fid,'%lg',[header.ncols header.nrows]);
fclose(fid);
Z(Z==header.nodata) = NaN;
Z = Z';
% create X and Y using meshgrid
refmat = [0 -header.cellsize;...
          header.cellsize 0;...
          header.xllcorner+(0.5*header.cellsize) - header.cellsize  ...
          (header.yllcorner+(0.5*header.cellsize))+((header.nrows)*header.cellsize)];

end



    
    
    

    
    