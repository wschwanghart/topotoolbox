function IX = getgriddedline(DEM,x,y)

%GETGRIDDEDLINE Rasterizes a line and returnes linear indices of the line
%
% Syntax 
%     
%     IX = getgriddedline(DEM,x,y)
%
% Description
%
%     getgriddedline takes a GRIDobj DEM, and two coordinate vectors as 
%     input, calculates the gridded path and returns a vector of linear 
%     indices IX into DEM.Z that describe the path. Note that the path must 
%     not intersect itself and should be not outside of the extent of the
%     DEM.
%
% Input arguments
%
%     DEM      GRIDobj
%     x,y      coordinate vectors
%
% Output arguments
%
%     IX       Linear index of the rasterized path.
%
%
% See also: line2GRIDobj, table2STREAMobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 07. June, 2023


%% find lines that are completely outside the grid boundaries
% force column vectors
x = x(:);
y = y(:);
xx = [x(1:end-1) x(2:end)];
yy = [y(1:end-1) y(2:end)];

lims = getextent(DEM);
I = all(xx > lims(2),2) | all(xx < lims(1),2) | ...
    all(yy > lims(4),2) | all(yy < lims(3),2);
I = [I;false];
I = ~I;
x = x(I);
y = y(I);

%%
[X,Y] = refmat2XY(DEM.refmat,DEM.size);
X = X(:);
Y = Y(:);

dx  = X(2)-X(1);
dy  = Y(2)-Y(1);

IX1 = (x-X(1))./dx + 1;
IX2 = (y-Y(1))./dy + 1;

cols = round(IX1);
rows = round(IX2);

%%
siz = DEM.size;
subs = [];

for r = 2:numel(rows)
    if any(isnan(rows([r r-1])))
        continue
    end
    
    p1 = [rows(r-1) cols(r-1)];
    p2 = [rows(r)   cols(r)  ];

    subsnew = getline(p1,p2);

    %if r == 2    
        subs = [subs;subsnew]; %#ok<AGROW>
    %else
    %    if size(subsnew,1) > 1
    %        subs = subsnew(2:end,:);
    %    end
    %end
end

if isempty(subs)
    return
end

% Remove cells that are outside the grid borders
outsidegrid = subs(:,1) < 1 | subs(:,1) > siz(1) | ...
              subs(:,2) < 1 | subs(:,2) > siz(2);
          
subs(outsidegrid,:) = [];

if isempty(subs)
    IX = [];
    return
end

IX   = sub2ind(siz,subs(:,1),subs(:,2));
IX   = unique(IX,'stable');

function subs = getline(p1,p2)

p = p1;
d = p2-p1;
N = max(abs(d));
s = d/N;

subs = zeros(N,2);
subs(1,:) = p1;

for ii=2:N+1
   p = p+s;
   subs(ii,:) = round(p);
end
end
end

