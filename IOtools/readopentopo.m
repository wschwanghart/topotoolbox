function DEM = readopentopo(varargin)

%READOPENTOPO read DEM using the opentopography.org API
%
% Syntax
%
%     DEM = readopentopo(pn,pv,...)
%
% Description
%
%     readopentopo reads DEMs from opentopography.org using the API
%     described on:
%     http://www.opentopography.org/developers
%     The DEM comes in geographic coordinates (WGS84) and should be
%     projected to a projected coordinate system (use reproject2utm) before
%     analysis in TopoToolbox.
%
%     NOTE: Starting on January 1st, 2022, an API authorization key will be 
%     required for this API. Users can request an API key via myOpenTopo in 
%     the OpenTopography portal (https://opentopography.org/developers).
%
% Input arguments
%
%     Parameter name values
%     'filename'       provide filename. By default, the function will save
%                      the DEM to a temporary file in the system's temporary 
%                      folder. The option 'deletefile' controls whether the
%                      file is kept on the hard drive.
%     'extent'         GRIDobj or four element vector with geographical 
%                      coordinates in the order [west east south north].
%                      If a GRIDobj is supplied, readopentopo uses the
%                      function GRIDobj/getextent to obtain the bounding
%                      box in geographical coordinates. If extent is set,
%                      then the following parameter names 'north',
%                      'south', ... are ignored.
%     'addmargin'      Expand the extent derived from 'extent',GRIDobj by a
%                      scalar value in °. Default is 0.01. The option is
%                      only applicable if extent is provided by a GRIDobj.
%     'north'          northern boundary in geographic coordinates (WGS84)
%     'south'          southern boundary
%     'west'           western boundary
%     'east'           eastern boundary
%     'demtype'        The global raster dataset *
%                      {'SRTMGL3'}:  SRTM GL3 (90m) (default)
%                      'SRTMGL1':    SRTM GL1 (30m)  
%                      'SRTMGL1_E':  SRTM GL1 (Ellipsoidal)  
%                      'AW3D30':     ALOS World 3D 30m  
%                      'AW3D30_E':   ALOS World 3D (Ellipsoidal)
%                      'SRTM15Plus': Global Bathymetry SRTM15+ V2.1 (only 
%                                    mediterranean area so far)
%                      'NASADEM':    NASADEM Global DEM 
%                      'COP30':      Copernicus Global DSM 30m 
%                      'COP90':      Copernicus Global DSM 90m 
%                        
%                      * requires API Key (see option 'apikey').
%
%     'apikey'         char or string. Users can request an API key via 
%                      myOpenTopo in the OpenTopography portal. You can
%                      also create a text file in the folder IOtools named
%                      opentopography.apikey which must contain the API
%                      Key. If there is a file, there is no need to provide
%                      this parameter name/value pair.
%     'verbose'        {true} or false. If true, then some information on
%                      the process is shown in the command window
%     'deletefile'     {true} or false. True, if file should be deleted
%                      after it was downloaded and added to the workspace.
% 
% Output arguments
%
%     DEM            Digital elevation model in geographic coordinates
%                    (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEM2 = readopentopo('extent',DEM);
%     DEM2 = reproject2utm(DEM2,90);
%     imagesc(DEM2)
%     hold on
%     getoutline(DEM)
%     hold off
%
% See also: GRIDobj, websave
%
% Reference: http://www.opentopography.org/developers
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 21. May, 2021


p = inputParser;
addParameter(p,'filename',[tempname '.tif']);
addParameter(p,'interactive',false);
addParameter(p,'extent',[]);
addParameter(p,'addmargin',0.01);
addParameter(p,'north',37.091337);
addParameter(p,'south',36.738884);
addParameter(p,'west',-120.168457);
addParameter(p,'east',-119.465576);
addParameter(p,'demtype','SRTMGL3');
addParameter(p,'deletefile',true);
addParameter(p,'verbose',true);
addParameter(p,'apikey',[]);
parse(p,varargin{:});

demtype = validatestring(p.Results.demtype,...
    {'SRTMGL3','SRTMGL1','SRTMGL1_E',...
     'AW3D30','AW3D30_E','SRTM15Plus',...
     'NASADEM','COP30','COP90'},'readopentopo');
%url = 'http://portal.opentopography.org/otr/getdem';
url = 'https://portal.opentopography.org/API/globaldem?';

% create output file
f = fullfile(p.Results.filename);

% check api

if isempty(p.Results.apikey)
    % check whether file opentopography.apikey is available
    if exist('opentopography.apikey','file')
        fid = fopen('opentopography.apikey');
        apikey = textscan(fid,'%c');
        apikey = apikey{1}';
        % Remove trailing blanks, if there are any
        apikey = deblank(apikey);
    else
        error('The DEM types NASADEM, COP30 and COP90 require an API Key')
    end
else
    apikey = p.Results.apikey;
end


% save to drive
options = weboptions('Timeout',100000);

% get extent
if ~isempty(p.Results.extent)
    if isa(p.Results.extent,'GRIDobj')
        ext = getextent(p.Results.extent,true);
        west  = ext(1) - p.Results.addmargin;
        east  = ext(2) + p.Results.addmargin;
        south = ext(3) - p.Results.addmargin;
        north = ext(4) + p.Results.addmargin;
        
    elseif numel(p.Results.extent) == 4
        west = p.Results.extent(1);
        east = p.Results.extent(2);
        south = p.Results.extent(3);
        north = p.Results.extent(4);
    else
        error('Unknown format of extent')
    end
else
    west = p.Results.west;
    east = p.Results.east;
    south = p.Results.south;
    north = p.Results.north;
end

% now we have an extent. Or did the user request interactively choosing
% the extent.
if any([isempty(west) isempty(east) isempty(south) isempty(north)]) || p.Results.interactive
    wm = webmap;
    % get dialog box
    messagetext = ['Zoom and resize the webmap window to choose DEM extent. ' ...
                         'Click the close button when you''re done.'];
    d = waitdialog(messagetext);
    uiwait(d);    
    [latlim,lonlim] = wmlimits(wm);
    west = lonlim(1);
    east = lonlim(2);
    south = latlim(1);
    north = latlim(2);
end
    
if p.Results.verbose
    a = areaint([south south north north],...
                [west east east west],almanac('earth','radius','kilometers'));
    disp('-------------------------------------')
    disp('readopentopo process:')
    disp(['DEM type: ' demtype])
    disp(['API url: ' url])
    disp(['Local file name: ' f])
    disp(['Area: ' num2str(a,2) ' sqkm'])
    disp('-------------------------------------')
    disp(['Starting download: ' datestr(now)])
end

% Download with websave
if isempty(apikey)
    outfile = websave(f,url,'west',west,...
        'east',east,...
        'north',north,...
        'south',south,...
        'outputFormat', 'GTiff', ...
        'demtype', demtype, ...
        options);
else
    outfile = websave(f,url,'west',west,...
        'east',east,...
        'north',north,...
        'south',south,...
        'outputFormat', 'GTiff', ...
        'demtype', demtype, ...
        'API_Key',apikey,...
        options);
end

if p.Results.verbose
    disp(['Download finished: ' datestr(now)])
    disp(['Reading DEM: ' datestr(now)])
end

try
    warning off
    DEM      = GRIDobj(f);
    warning on
    
    msg = lastwarn;
    if ~isempty(msg)
        disp(' ')
        disp(msg)
        disp(' ')
    end
    
    DEM.name = demtype;
    if p.Results.verbose
        disp(['DEM read: ' datestr(now)])
    end
    
catch
    % Something went wrong. See whether we can derive some information.
    fid = fopen(outfile);
    in = textscan(fid,'%c');
    disp('Could not retrieve DEM. This is the message returned by opentopography API:')
    disp([in{1}]')
    disp('readopentopo returns empty array')
    fclose(fid);
    DEM = [];
end
    

if p.Results.deletefile
    delete(f);
    if p.Results.verbose
        disp('Temporary file deleted')
    end
end

if p.Results.verbose
    disp(['Done: ' datestr(now)])
    disp('-------------------------------------')
end
end

function d = waitdialog(messagetext)
    d = dialog('Position',[300 300 250 150],'Name','Choose rectangle region',...
        'WindowStyle','normal');

    txt = uicontrol('Parent',d,...
               'Style','text',...
               'Position',[20 80 210 40],...
               'String',messagetext);

    btn = uicontrol('Parent',d,...
               'Position',[85 20 70 25],...
               'String','Close',...
               'Callback','delete(gcf)');
end


