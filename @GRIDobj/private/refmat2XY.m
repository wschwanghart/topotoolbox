function [x,y] = refmat2XY(R,siz)

% Convert referencing matrix (refmat) to coordinate vectors
%
% Syntax
%
%     [x,y] = refmat2XY(R,siz)
%
% Input arguments
%
%     R     referencing matrix
%     siz   size of the grid
%
% Output arguments
%
%     x     x-coordinates
%     y     y-coordinates
%
%

nrrows = siz(1);
nrcols = siz(2);

x = [ones(nrcols,1) (1:nrcols)' ones(nrcols,1)]*R;
x = x(:,1)';

y = [(1:nrrows)' ones(nrrows,2)]*R;
y = y(:,2);

