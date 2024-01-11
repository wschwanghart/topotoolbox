function [S] = stat(SW,varargin)

%STAT calculates statistics for swath profile object (SWATHobj)
%
% Syntax
%
%    S = stat(SW)
%    S = stat(SW,type)
%
% Description
%
%     stat(SW) calculates the arithmetic mean across a SWATHobj for each
%     point along it.
%
%     stat(SW,type) calculates a statistical metric across a SWATHobj for 
%     each point along it. type specifies the metric, with valid strings
%     being {'min','max','mean','range','prctile'}. When 'prctile' is
%     chosen, a cell with {'percentile',value} should be provided, with the
%     value being equal to the percentile.
%
% Input arguments
%
%     SW     swath object (Class: SWATHobj)
%     type   string specifying statistical metric 
%            {'min','max','mean' (default),'range','prctile'}
%
%
% Output
%
%     S      vector with statistic metric values along the SWATHobj
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: 2017
    


narginchk(1,2)

if ~isa(SW,'SWATHobj')
    error('TopoToolbox: First input argument must be of calss SWATHobj')
end

validtypes = {'min','max','mean','range','prctile'};
if nargin>1
    type = varargin{1};
    if iscell(type)
        value = type{2};
        type = type{1};
    end
    validatestring(type,validtypes);
else
    type = 'mean';
end


switch type
    case 'prctile'
        S = prctile(SW.Z,value);
    otherwise
        str = sprintf('%s(SW.Z);',type);
        S = eval(str);
end





end

