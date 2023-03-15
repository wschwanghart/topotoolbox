function [S] = stat(SW,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


narginchk(1,2)

if ~isa(SW,'SWATHobj')
    error('TopoToolbox: First input argument must be of calss SWATHobj')
end

validtypes = {'min','max','mean','range','prctile','all'};
if nargin>1
    type = varargin{1};
    if iscell(type)
        value = type{2};
        type = type{1};
    end
    validatestring(type,validtypes);
end


switch type
    case 'prctile'
        S = prctile(SW.Z,value);
    otherwise
        str = sprintf('%s(SW.Z);',type);
        S = eval(str);
end





end

