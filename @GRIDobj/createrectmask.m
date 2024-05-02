function MASK = createrectmask(DEM,usehillshade)

%CREATEMASK create a binary mask using rectangle mapping
%
% Syntax
%
%     MASK = createrectmask(DEM)
%     MASK = createrectmask(DEM,usehillshade)
%
% Description
%
%     createrectmask is an interactive tool to create a mask based on an
%     interactively mapped axis-aligned rectangle.
%
% Input arguments
%
%     DEM           GRIDobj
%     usehillshade  use hillshade as background image ({false} or true)
%
% Output arguments
%
%     MASK   GRIDobj with logical mask
%
%
% See also: imroi, GRIDobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 2. May, 2024

if nargin == 1
    usehillshade = false;
end

figure
if usehillshade
   imageschs(DEM);
else
   imagesc(DEM);
end
title('Create polygon');
usedrawpolygon = ~verLessThan('matlab','9.5');
if usedrawpolygon
    ext = getextent(DEM);  

    h = drawrectangle('DrawingArea',...
        [ext(1) ext(3) ext(2)-ext(1) ext(4)-ext(3)]);
    pos = customWait(h);
else
    h = imrect; 
    pos = wait(h);
end

MASK = DEM;
MASK.name = 'mask';
MASK.Z = createMask(h);

if usedrawpolygon
    pos = h.Position;
else
    pos = getPosition(h);
end
delete(h);
hold on
plot(pos([1:end 1],1),pos([1:end 1],2));
hold off
title('done');


end


function pos = customWait(hROI)

% Listen for mouse clicks on the ROI
l = addlistener(hROI,'ROIClicked',@clickCallback);

% Block program execution
uiwait;

% Remove listener
delete(l);

% Return the current position
pos = hROI.Position;

end

function clickCallback(~,evt)

if strcmp(evt.SelectionType,'double')
    uiresume;
end

end
