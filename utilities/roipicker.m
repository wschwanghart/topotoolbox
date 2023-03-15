function ext = roipicker(varargin)

%ROIPICKER Choose region in geoaxes
%
% Syntax
%
%     ext = roipicker
%     ext = roipicker(pn,pv)
%
% Description
%
%     roipicker opens a GUI which lets user choose a rectangle extent in
%     a global map. The function stops execution until an extent is
%     retrieved and returns the extent. 
%
%     If users quit selecting a region by closing one of the windows, ext
%     will be empty.
%
% Input arguments
%
%     Following parameter name/value combinations apply
%
%     'basemap'        Basemap used in map. Default is 'topographic'. See
%                      geoaxes for alternatives. Other basemaps can be
%                      selected within the roipicker GUI.
%     'ellipsoid'      ellipsoid used for calculating the area of the ROI.
%                      Default is referenceSphere('earth','km'). There will
%                      rarely be a need to change this parameter.
%     'ask'            {true} or false. If true, users will have to confirm 
%                      if they quit the application without selecting a ROI.
%     'drawingarea'    By default, the drawing area of the ROI is limited to 
%                      [-90 -180 180 360]. Adjust the extent with a vector
%                      containing [latmin lonmin height width].
%     'requestlimit'   {inf}. Provide a maximum area in km here. If the ROI
%                      exceeds the request limit, the select button is
%                      deactivated and the area label is red.
%
% Output arguments
%
%     ext   [minlat minlon maxlat maxlon]
%
% See also: readopentopo, readopenalti
%
% Authors:  Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
%
% Date: 28. September, 2022

% Parse inputs 
p = inputParser;
p.FunctionName = 'roipicker';
addParameter(p,'basemap','topographic')
addParameter(p,'ellipsoid',referenceSphere('earth','km'));
addParameter(p,'ask',true)
addParameter(p,'drawingarea',[-90 -180 180 360])
addParameter(p,'requestlimit',inf) % in km
parse(p,varargin{:})

% Define ext 
ext = [];

% Control figure
cfig = uifigure('Units','Normalized','Position',[0 .60 .20 .30],'Name','Controls',...
    'NumberTitle','off');
% Map figure
mfig = figure('Units','Normalized','Position',[0.20 .1 .70 .80],...
    'ToolBar','none','MenuBar','none','NumberTitle','off',...
    'Name','Map');

% Create map
warning off
gx = geoaxes(mfig,'Basemap',p.Results.basemap);
gx.Toolbar.Visible = 'on';
enableDefaultInteractivity(gx)

% Roi is an empty array at the beginning but will be a roi-object later
% which spans several functions
roi = [];

% Grid layout for controls
grid1 = uigridlayout(cfig,[8 1]);
grid1.RowHeight = {22,22,22, 22, 22, 22};
grid1.ColumnWidth = {'1x'};

% Create button for rectangle selection
b1 = uibutton(grid1,'state','Text','Draw ROI',...
    'ValueChangedFcn', @(btn,event) drawrect(btn,gx),...
    'Tooltip','Start drawing a rectangle.');
b2 = uibutton(grid1,'Text','Clear ROI',...
    'ButtonPushedFcn', @(btn,event) clearroi(btn,gx),...
    'Tooltip','Delete current ROI.');
b3 = uibutton(grid1,'Text','Set ROI to map extent',...
    'ButtonPushedFcn', @(btn,event) mapextent2roi(btn,gx),...
    'Tooltip','Set ROI to the current extent of the map.');
b4 = uibutton(grid1,'Text','Zoom to ROI',...
    'ButtonPushedFcn', @(btn,event) zoom2roi(btn,gx),...
    'Tooltip','Zoom map to the extent of the ROI.',...
    'Enable',true);
b5 = uibutton(grid1,'text','Select and close',...
    'ButtonPushedFcn', @(btn,event) select(btn,gx),'UserData',0,...
    'Tooltip','Retrieve current extent of ROI and close figures.',...
    'Enable', false);

% Dropdown list for choosing basemaps
items = {'topographic','satellite','landcover','streets-light',...
         'streets-dark','streets','grayterrain',...
         'colorterrain', 'bluegreen', 'darkwater','grayland'};

dd1 = uidropdown(grid1,'Editable','off','Items',items,...
    'ValueChangedFcn',@(dd,event) basemapselect(dd,p));

% Labels to display information about
lbl1 = uilabel(grid1,"Text",'Area = NaN', 'Interpreter','html',...
    'FontName','Arial');
lbl2 = uilabel(grid1,"Text",'Extent = NaN');

if p.Results.ask
cfig.CloseRequestFcn = @my_closereq;
mfig.CloseRequestFcn = @my_closereq;
end

waitfor(b5,'UserData');



% ----------------------------
    function drawrect(btn,gx)
        % Callback: Draw rectangle

        disableDefaultInteractivity(gx)

        if btn.Value
            % check if roi already exists
            if ~isempty(roi)
                delete(roi)
                roi = [];
            end

            % Deactivate all buttons
            b1.Enable = false;
            b2.Enable = false;
            b3.Enable = false;
            b4.Enable = false;
            b5.Enable = false;

            roi = drawrectangle(gx,'DrawingArea',p.Results.drawingarea,...
                'Deletable', false);
            displayarea
            addlistener(roi,'MovingROI',@allevents);
            addlistener(roi,'ROIMoved',@allevents);

            % Activate all buttons
            b1.Enable = true;
            b2.Enable = true;
            b3.Enable = true;
            b4.Enable = true;
            b5.Enable = true;

        else
            if exist('roi','var')
                delete(roi)
                roi = [];
            end
        end
        enableDefaultInteractivity(gx);
    end
    
    % --------------
    function clearroi(btn,gx)
        % Callback: Delete ROI

        b1.Value = 0;
        if isempty(roi)
            return
        else
            delete(roi)
            roi = [];
        end
        b5.Enable = false;
        b4.Enable = false;
    end

    % --------------
    function zoom2roi(btn,gx)
        % Callback: Zoom to ROI
        if isempty(roi)
            return
        end

        % get position of roi
        pos = roi.Position;
        pos(3) = pos(3)+pos(1);
        pos(4) = pos(4)+pos(2);

        geolimits(gx,pos([1 3]),pos([2 4]))
    end

    % --------------
    function mapextent2roi(btn,gx)
        % Callback: Adjust mapextent to match ROI

        % Delete an existing roi
        if exist('roi','var')
            delete(roi)
        end
        % Get mapextent
        [latLimits,lonLimits] = geolimits(gx);

        lonLimits(1) = min(max(lonLimits(1),-180), 180);
        lonLimits(2) = min(max(lonLimits(2),-180), 180);

        if lonLimits(1) == lonLimits(2)
            return
        end

        roi = drawrectangle(gx,'Position',[latLimits(1) lonLimits(1) ...
            latLimits(2)-latLimits(1) lonLimits(2)-lonLimits(1)],...
            'DrawingArea',p.Results.drawingarea);
        displayarea
        addlistener(roi,'MovingROI',@allevents);
        addlistener(roi,'ROIMoved',@allevents);

        b4.Enable = true;
        b5.Enable = true;
    end

    % --------------
    function allevents(src,evt)
        % Callback: Adjust behavior of geoaxes when ROI is moved
        evname = evt.EventName;

        switch(evname)
            case{'MovingROI'}
                disableDefaultInteractivity(gx)
                displayarea
                
            case{'ROIMoved'}
                enableDefaultInteractivity(gx)
                displayarea
                
        end
    end

    % --------------
    function select(src,evt)
        % Callback: Select and output ROI extent and finish
        if isempty(roi)
            return
        end
        pos = roi.Position;
        pos(3) = pos(3)+pos(1);
        pos(4) = pos(4)+pos(2);
        ext = pos;

        b5.UserData = 1;

        cfig.CloseRequestFcn = 'closereq';
        mfig.CloseRequestFcn = 'closereq';

        close(cfig)
        close(mfig)

        warning on
    end

    % --------------
    function basemapselect(dd,evt)
        % Callback: Select basemap
        gx.Basemap = dd.Value;

    end

    % --------------
    function displayarea
        % Callback: Display area and extent of ROI
        pos = roi.Position;
        pos(3) = pos(3)+pos(1);
        pos(4) = pos(4)+pos(2);

        area    = areaquad(pos(1),pos(2),pos(3),pos(4),p.Results.ellipsoid);
        areastr = ['Area: ' num2str(area,2) ' km<sup>2</sup>'];
        
        if area > p.Results.requestlimit
            lbl1.FontColor = 'r';
            b5.Enable = false;
            areastr = [areastr ' (ROI too large)'];
        else
            lbl1.FontColor = 'k';
            b5.Enable = true;
        end
        lbl1.Text = areastr;

        posstr = ['T:' num2str(pos(3)) ' B:' num2str(pos(1)) ...
            ' L:' num2str(pos(2)) ' R:' num2str(pos(4))];
        lbl2.Text = posstr;

    end

    % --------------
    function my_closereq(src,event)        
        % Close request function
        % to display a question dialog box
        selection = questdlg('Quit selection and close?', ...
            'Close Request Function', ...
            'Yes','No','Yes');
        switch selection
            case 'Yes'
                delete(cfig)
                delete(mfig)
            case 'No'
                return
        end
    end

end
