function tf = validatealignment(GRID1,GRID2)

%VALIDATEALIGNMENT validates whether instances of GRIDobj are spatially aligned
%
% Syntax
%
%     tf = validatealignment(GRID1,GRID2)
%     validatealignment(GRID1,GRID2)
%
% Description
%
%     returns true if instances of GRIDobj are spatially 
%     aligned. When the function returns false and is called without 
%     output argument, the function returns an error message.
%
% Input arguments
%
%     GRID1 instance of GRIDobj
%     GRID2 instance of GRIDobj or matrix
%
% Output arguments
% 
%     tf    true or false
%
% See also: GRIDobj
%   
% Author:  Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 8. August, 2015 

% check if geometric properties of a FLOWobj and GRIDobj instance are equal
if isa(GRID2,'GRIDobj')
    TF = isequal(GRID1.size,GRID2.size) && isequal(GRID1.refmat,GRID2.refmat);
else
    TF = isequal(GRID1.size,size(GRID2));
end

if nargout == 1;
    tf = TF;
else
    if ~TF
        if isa(GRID2,'GRIDobj')
            error('TopoToolbox:incorrectinput',...
                ['The two GRIDobj instances do not align each other. Make sure that \n' ...
                'both instances have the same spatial reference. Both variables \n' ...
                'are deemed to have the same reference if their properties ''size'' \n' ...
                'and ''refmat'' are both equal.']);
        else
            error('TopoToolbox:incorrectinput',...
                ['GRIDobj and input matrix do not align each other. Make sure that \n' ...
                'both instances have the same spatial reference. Both variables \n' ...
                'are deemed to have the same reference if the GRIDobj''s property \n' ...
                '''size'' and size(A) is equal.']);
        end
    end
end
