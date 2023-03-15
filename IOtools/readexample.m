function DEM = readexample(example,varargin)

%READEXAMPLE read DEM from TopoToolbox DEM github repository
%
% Syntax
%
%     DEM = readexample(example)
%
% Description
%
%     readexample reads DEMs from the TopoToolbox DEM repository.
%     Currently available examples are
%
%     'kunashiri'
%     'perfectworld'
%     'taalvolcano'
%     'taiwan'
%     'tibet'
%
%     The function requires internet connection.
%
% Input arguments
%
%     example    string of example name (e.g. 'tibet')
%
% Output arguments
%
%     DEM        example data. This is usually a GRIDobj. However, in some
%                cases, DEM is a structure array with variable data.
%     
% See also: GRIDobj, websave, readopentopo
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. March, 2020



example = lower(example);

p = inputParser;
addParameter(p,'filename',[tempname '.tif']);
addParameter(p,'deletefile',true);
addParameter(p,'verbose',true);
parse(p,varargin{:});

% create output file
f = fullfile(p.Results.filename);

switch example
    case 'taiwan'
        url = 'https://github.com/wschwanghart/DEMs/raw/master/taiwan.tif';
        istif = true;
    case 'tibet'
        url = 'https://github.com/wschwanghart/DEMs/raw/master/tibet.tif';
        istif = true;
    case 'taalvolcano'
        url = 'https://github.com/wschwanghart/DEMs/raw/master/taalvolcano.tif';
        istif = true;
    case 'kunashiri'
        url = 'https://github.com/wschwanghart/DEMs/raw/master/kunashiri.tif';
        istif = true;
    case 'perfectworld'
        url = 'https://github.com/wschwanghart/DEMs/raw/master/perfectworld.tif';
        istif = true;
    otherwise 
        error('There is no such example file.')
end
      
    

% save to drive
options = weboptions('Timeout',100000);

% Download with websave
if p.Results.verbose
    disp([datestr(now) ' -- Downloading...'])
end
outfile = websave(f,url,options);
if p.Results.verbose
    disp([datestr(now) ' -- Download finished...'])
end

% Read grid or mat file
if istif
    DEM = GRIDobj(f);
    DEM.name = example;
else
    DEM = load(f); 
end

if p.Results.deletefile
    delete(f);
    if p.Results.verbose
        disp([datestr(now) ' -- Temporary file deleted'])
    end
end