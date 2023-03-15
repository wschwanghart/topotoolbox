function tf = validatealignment(FD,G)

%VALIDATEALIGNMENT validates whether instances of FLOWobj and GRIDobj are spatially aligned
%
% Syntax
%
%     tf = validatealignment(FD,G)
%     validatealignment(FD,G)
%
% Description
%
%     returns true if instances of FLOWobj and GRIDobj are spatially 
%     aligned. When the function returns false and is called without 
%     output argument, the function returns an error message.
%
% Input arguments
%
%     FD    instance of FLOWobj
%     G     instance of GRIDobj
%
% Output arguments
% 
%     tf    true or false
%
% See also: FLOWobj, GRIDobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. March, 2016

% check if geometric properties of a FLOWobj and GRIDobj instance are equal
if isa(G,'GRIDobj')
    TF = isequal(FD.size,G.size) && isequal(FD.refmat,G.refmat);
else
    TF = isequal(FD.size,size(G));
end

if nargout == 1
    tf = TF;
else
    if ~TF
        if isa(G,'GRIDobj')
            error('TopoToolbox:incorrectinput',...
                ['FLOWobj and GRIDobj do not align each other. Make sure that \n' ...
                'both instances have the same spatial reference. Both variables \n' ...
                'are deemed to have the same reference if their properties ''size'' \n' ...
                'and ''refmat'' are both equal.']);
        else
            error('TopoToolbox:incorrectinput',...
                ['FLOWobj and input matrix do not align each other. Make sure that \n' ...
                'both instances have the same spatial reference. Both variables \n' ...
                'are deemed to have the same reference if the FLOWobj''s property \n' ...
                '''size'' and size(A) is equal.']);
        end
    end
end
