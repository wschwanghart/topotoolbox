function G = STREAMobj2GRIDobj(S,nal)

%STREAMOBJ2GRIDOBJ convert an instance of STREAMobj to an instance of GRIDobj
%
% Syntax 
%
%     G = STREAMobj2GRIDobj(S)
%     G = STREAMobj2GRIDobj(S,nal)
%
% Description
%     
%     STREAMobj2GRIDobj converts an instance of STREAMobj S to a new instance
%     of GRIDobj G with the same spatial reference and extent as the instance
%     from which S has been derived. If the second input is a
%     node-attribute list, STREAMobj2GRIDobj will write the node values to
%     the respective locations in the output GRIDobj G.
%     
% Input arguments
%
%     S     STREAMobj
%     nal   node attribute list
%
% Output arguments
%
%     G     GRIDobj (contains logical raster with ones at stream locations
%           and zeros elsewhere. If a node attribute list was supplied, G
%           will have the same underlying class as nal. All non-stream
%           values are set to nan.
%
%
% See also: STREAMobj, GRIDobj
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. May, 2016

narginchk(1,2)

if nargin == 1;
    nal = true(size(S.x));
else
    if ~isnal(S,nal);
        error('TopoToolbox:isnal',...
            'The second input argument is not a valid node attribute list');
    end
end

G = GRIDobj([]);

props_G = fieldnames(G);
props_S = fieldnames(S);

[props_joint,b] = ismember(props_S,props_G); 

six = find(props_joint);
gix = b(six);

for r = 1:numel(six);
    G.(props_G{gix(r)}) = S.(props_S{six(r)});
end


if islogical(nal)
    G.Z = false(G.size);
elseif isinteger(nal)
    G.Z = zeros(G.size,'like',nal);
elseif isfloat(nal); 
    G.Z = nan(G.size,'like',nal);
end
G.Z(S.IXgrid) = nal;

    

G.name = 'stream grid';





