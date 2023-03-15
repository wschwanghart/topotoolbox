function GRIDobj2geotiff(A,file)

%GRIDobj2geotiff Exports an instance of GRIDobj to a geotiff file
%
% Syntax
%    
%     GRIDobj2geotiff(DEM)
%     GRIDobj2geotiff(DEM,filename)
%
% Description
%
%     GeoTIFF is a common image file format that stores coordinates and
%     projection information to be read by most GIS software.
%     GRIDobj2geotiff writes an instance of GRIDobj to a GeoTIFF file. 
%
%     GRIDobj2geotiff requires the function geotiffwrite available with 
%     the Mapping Toolbox. If geotiffwrite does not exist on the search
%     path, the function will write a standard tif together with a
%     '.tfw'-file (worldfile, http://en.wikipedia.org/wiki/World_file ) to
%     the disk. 
%
%     GRIDobj2geotiff(DEM) opens a dialogue box to save the GeoTIFF
%
%     GRIDobj2geotiff(DEM,filename) saves the DEM to the specified 
%     filename
%
% Input arguments
%
%     DEM        instance of GRIDobj
%     filename   absolute or relative path and filename
%
% See also: GRIDobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017

narginchk(1,2)

% if only 1 argument, open file dialog box
if nargin == 1
    [FileName,PathName] = uiputfile({'*.tif'});
    if FileName == 0
        disp('     no output written to disk')
        return
    end
    file = [PathName FileName];
end

% try to use geotiffwrite, which comes with the Mapping Toolbox
try
    if isempty(A.georef);
        geotiffwrite(file,A.Z,A.refmat);
    else
        geotiffwrite(file,A.Z,A.georef.SpatialRef,...
            'GeoKeyDirectoryTag',A.georef.GeoKeyDirectoryTag);
    end
catch ME
    warning('TopoToolbox:GRIDobj',...
        ['GRIDobj2geotiff is unable to write a geotiff. Either you don''t \n'...
         'have the mapping toolbox, or there was another issue with geotiffwrite. \n'...
         'GRIDobj2geotiff instead writes a tif-image together with a world \n'...
         'file (*.tfw) which contains data on spatial referencing of the \n' ...
         'image, yet which lacks information on the type of projection used.']); 
         
         
    % if geotiffwrite is not available or any other error occurs
    % a tif file will be written to the disk together with a worldfile
    % .tfw-file.
    [pathstr, name, ~] = fileparts(file);
    k = refmat2worldfile(A.refmat);
    dlmwrite(fullfile(pathstr,[name '.tfw']),k,'precision', '%.10f');
    A = A.Z;
    
    
    siz = size(A);
    cla = class(A);
    
    switch cla;
        case 'double'
            BpS = 64;
            TSF = Tiff.SampleFormat.IEEEFP;
        case 'single'
            BpS = 32;
            TSF = Tiff.SampleFormat.IEEEFP;
        otherwise
            if islogical(A);
                A = uint32(A);
                cla = 'uint32';
            end
            BpS = round(log2(double(intmax(cla))));
            TSF = Tiff.SampleFormat.UInt;
            
    end
    
    t = Tiff(file,'w');
    tagstruct.ImageLength = siz(1);
    tagstruct.ImageWidth = siz(2);
    tagstruct.BitsPerSample = BpS;
    tagstruct.SampleFormat = TSF;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    tagstruct.Photometric = 0;
    
    t.setTag(tagstruct);
    t.write(A);
    t.close;
end

end

function k = refmat2worldfile(r)
% does not support rotation

k(1,1) = r(2,1);
k(4,1) = r(1,2);
k(5,1) = r(3,1)+k(1);
k(6,1) = r(3,2)+k(4);
end