function tf = validatealignment(S,G)

%VALIDATEALIGNMENT is an instance of STREAMobj is spatially aligned with another object of TopoToolbox
%
% Syntax
%
%     tf = validatealignment(S,G)
%     validatealignment(S,G)
%
% Description
%
%     returns true if instances of STREAMobj and GRIDobj are spatially 
%     aligned. When the function returns false and is called without 
%     output argument, the function returns an error message.
%
% Input arguments
%
%     S     instance of STREAMobj
%     G     instance of GRIDobj, STREAMobj or FLOWobj
%
% Output arguments
% 
%     tf    true or false
%
%
% See also: STREAMobj
%   
% Author:  Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 8. August, 2015 

% check if geometric properties of a FLOWobj and GRIDobj instance are equal
if isa(G,'GRIDobj') || isa(G,'STREAMobj') || isa(G,'FLOWobj');
    TF = isequal(S.size,G.size) && isequal(S.refmat,G.refmat);
else
    TF = isequal(S.size,size(G));
end

if nargout == 1;
    tf = TF;
else
    if ~TF
        if isa(G,'GRIDobj')
            error('TopoToolbox:incorrectinput',...
                ['STREAMobj and GRIDobj do not align each other. Make sure that \n' ...
                'both instances have the same spatial reference. Both variables \n' ...
                'are deemed to have the same reference if their properties ''size'' \n' ...
                'and ''refmat'' are both equal.']);
        else
            error('TopoToolbox:incorrectinput',...
                ['STREAMobj and input matrix do not align each other. Make sure that \n' ...
                'both instances have the same spatial reference. Both variables \n' ...
                'are deemed to have the same reference if the FLOWobj''s property \n' ...
                '''size'' and size(A) is equal.']);
        end
    end
end
