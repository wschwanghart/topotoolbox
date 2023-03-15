function h = plotdz(SW,varargin)
%PLOTDZ creates distance-elevation plot of SWATHobj
%
% Syntax
%
%     plotdz(SW)
%     plotdz(SW,'pn','pv',...)
%     h = ...;
%
% Description
%
%     PLOTDZ creates a profile-view plot of a SWATHobj, showing the
%     statictics (min, max, mean +/- standard dev.) of the z-values as a
%     function of distance along the profile of the SWATHobj.
%
% Input arguments
%
%     SW     instance of SWATHobj
%
%     Parameter name/value pairs
%
%     'left'   {true}, false
%     determines if the left half (as seen from the direction of the 
%     central line) of the SWATHobj is included in the statistics
%
%     'right'   {true}, false
%     determines if the right half (as seen from the direction of the 
%     central line) of the SWATHobj is included in the statistics
%
%     'distadjust'   {0}, scalar
%     allows shifting the x-axis by a scalar value. This is useful when
%     alligning a SWATHobj that was obtained along, e.g., a reach of
%     a drainage network, with the STREAMobj that it was derived from
%
%     'boundedline'   true, {false}
%     allows the plot to be created with the 'boundedline' plotting
%     function by Kelley Kearny, if available (downloadable from Matlab 
%     Central website)
%
%
% Output arguments (optional)
%
%     h    handle object of the plotted features. h contains several
%          handles to the different parts of the plot
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     SW = SWATHobj(DEM,'dx',200,'dy',200);
%     G = gradient8(DEM,'degree');
%     SWG = mapswath(SW,G);
%     figure, plotdz(SWG)
%     title('Slope angle along swath')
%
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: June, 2013


% Parse inputs
p = inputParser;
p.FunctionName = 'plotdz';
addRequired(p,'SW',@(x) isa(x,'SWATHobj'));
addParamValue(p,'left',true,@(x) islogical(x))
addParamValue(p,'right',true,@(x) islogical(x))
addParamValue(p,'distadjust',0,@(x) isnumeric(x))
addParamValue(p,'boundedline',false,@(x) islogical(x))
parse(p,SW,varargin{:});

% parameters
distadjust = p.Results.distadjust;
left       = p.Results.left;
right      = p.Results.right;
boundedl   = p.Results.boundedline;

if ~left && ~right
    warning('Empty SWATHobj: nothing to plot')
    return;
end


if ~left
    ny = ceil(length(SW.disty)/2)+1;
    SW.Z = SW.Z(ny:end,:);
elseif ~right
    ny = floor(length(SW.disty)/2);
    SW.Z = SW.Z(1:ny,:);
end

z_min = nanmin(SW.Z,[],1);
z_max = nanmax(SW.Z,[],1);
z_mean = nanmean(SW.Z,1)';
z_std = nanstd(SW.Z,0,1)';
dist = SW.distx+distadjust;

if exist('boundedline','file') && (boundedl)
    % Use plotting function 'boundedline', by Kelley Kearny, if
    % found in search path. Available from Matlab Central.
    [hp(1), hp(2)] = boundedline(dist,z_mean,z_std,'alpha'); hold on
    hp(3) = plot(dist,z_min,'c-');
    plot(dist,z_max,'c-')
else
    hp(1) = plot(dist,z_mean,'r-'); hold on
    hp(2) = plot([dist;nan;dist],[z_mean-z_std;nan;z_mean+z_std],'b-');
    hp(3) = plot([dist;nan;dist],[z_min nan z_max],'c-');
end

drawnow
xlabel(sprintf('Distance along profile (%s)',SW.xyunit))
ylabel(sprintf('Z (%s)',SW.zunit))
legend(hp,{'Mean','+/- St.Dev.','Min/Max'})

if nargout == 1;
    h = hp;
end

