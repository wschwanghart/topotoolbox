function DEM = reclassify(DEM,varargin)

%RECLASSIFY generate univariate class intervals for an instance of GRIDobj
%
% Syntax
%
%     C = reclassify(DEM);
%     C = reclassify(DEM,'method',value)
%
% Description
%
%     reclassify bins continous values of an instance of GRIDobj by setting
%     class intervals based on different classification methods. The
%     default method is equal interval classification with 10 classes.
%
% Input Arguments
%
%     DEM       Grid (class = GRIDobj)
%               
%     'Methods' and [values]
%
%     'equalintervals'   [number of classes]
%     'definedintervals' [vector of class breaks]
%     'equalquantiles'   [number of classes]
%     'definedquantiles' [vector of quantiles] 
%                        e.g. [0.1 0.9] results in three classes                        
%     'kmeans'           [number of classes]
%                        uses the k_means algorithm of Yi Cao
%                        FEX submission 19344 included in this function
%                        ! may take a while for large grids !
%     'std'              [s = scale of standard deviation, e.g. 2]
%                        equal intervals with a width of std/s. Bins are
%                        centered around sample mean
%     'otsu'             [number of classes, must be (2 4 8 16 etc)] 
%                        recursive Otsu thresholding (see function 
%                        graythresh)
%
% Output argument
%
%     C        Classified grid (class = GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     C   = reclassify(DEM,'equalquantiles',10);
%     imageschs(DEM,C)
%
% See also: k_means, graythresh
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. January, 2013

narginchk(1,3)

allowedmethods = {'equalintervals',...
                  'definedintervals',...
                  'definedquantiles',...
                  'equalquantiles',...
                  'kmeans',...
                  'standarddeviation',...
                  'std',...
                  'otsu'...
                 };
% nr of classes
if nargin == 1;
    method = 'equalinterval';
    num = 10;
else
    method = validatestring(varargin{1},allowedmethods);
    
    if nargin==3
        num = varargin{2};
    end
end

INAN = isnan(DEM.Z) | isinf(DEM.Z);
DEM.Z(INAN) = nan;
inan = any(INAN(:));

switch method
    case 'equalintervals'
        validateattributes(num,{'numeric'},{'scalar'});
        
        if inan
            DEM.Z(~INAN) = mat2gray(DEM.Z(~INAN));
            DEM.Z = grayslice(DEM.Z,num);
            if num == 256;
                DEM.Z = double(DEM.Z) + 1;
                DEM.Z(INAN) = nan;
            elseif num < 256;
                DEM.Z = DEM.Z + 1;
                DEM.Z(INAN) = 0;
            else
                DEM.Z(INAN) = nan;
            end
            
        else
            DEM.Z = mat2gray(DEM.Z);
            DEM.Z = grayslice(DEM.Z,num);
            
            if num <= 255
                DEM.Z = DEM.Z + 1;
            elseif num == 256
                DEM.Z = double(DEM.Z) + 1;
            end
        end
    case 'definedintervals'
        validateattributes(num,{'numeric'},{'vector'});
        num(end+1) = inf;
        num = [-inf; num(:)];
        siz = size(DEM.Z);
        [~,DEM.Z] = histc(DEM.Z(:),num);
        DEM.Z = reshape(DEM.Z,siz);
    case 'equalquantiles'
        validateattributes(num,{'numeric'},{'scalar','integer','>',1});
        z = DEM.Z(~INAN);
        z = z(:);
        z = sort(z,'ascend');
        q = (1:num)/num;
        q = ceil(q*numel(z));
        edges = z(q);
        edges(end) = inf;
        [~,DEM.Z] = histc(DEM.Z(:),[-inf; edges]);
        DEM.Z = reshape(DEM.Z,DEM.size);
        DEM.Z(INAN) = nan;
    case 'definedquantiles'
        validateattributes(num,{'numeric'},{'vector','>',0,'<=',1});
        p = num(:);
        if p(end) < 1;
            p(end+1) = 1;
        end
        
        z = DEM.Z(~INAN);
        z = z(:);
        n = numel(z);
        edges = interp1((0:n-1).'/(n-1), sort(z), p);
  
        edges(end) = inf;
        [~,DEM.Z] = histc(DEM.Z(:),[-inf; edges]);
        DEM.Z = reshape(DEM.Z,DEM.size);
        DEM.Z(INAN) = nan;
        
    case {'standarddeviation','std'}
        validateattributes(num,{'numeric'},{'vector','>',0,'<=',1});
        z = DEM.Z(~INAN);
        z = z(:);
        s = std(z(:));
        s = s*num;
        mz = mean(z);
        minz = min(z);
        maxz = max(z);
        
        edges1 = [mz-s/2 :-s: minz];
        if edges1(end) > minz
            edges1(end+1) = -inf;
        else
            edges1(end) = -inf;
        end
        edges2 = mz+s/2 :s: maxz;
        
        edges2(end) = inf;
        
        [~,DEM.Z] = histc(DEM.Z(:),[edges1(end:-1:1) edges2]);
        DEM.Z = reshape(DEM.Z,DEM.size);
        DEM.Z(INAN) = nan;
        
    case 'kmeans';
        
        validateattributes(num,{'numeric'},{'scalar','integer','>',1});
        z = DEM.Z(~INAN);
        z = z(:);
        
        IX = k_means(z,num);
        
        DEM.Z(:,:) = nan;
        DEM.Z(~INAN) = IX;
    case 'otsu';        
        validateattributes(num,{'numeric'},{'scalar','integer','>',1});
        if ceil(log2(num)) ~= log2(num);
            error('TopToolbox:GRIDobj','log2(value) must be an integer')
        end
                
        INAN = ~INAN;
        z  = mat2gray(DEM.Z(INAN));
        DEM.Z(INAN) = cast(otsu(z,num),class(DEM.Z));
end

end




% subfunctions

function [gIdx,c]=k_means(X,k)
% K_MEANS    k-means clustring
%   IDX = k_means(X, K) partititions the N x P data matrix X into K
%   clusters through a fully vectorized algorithm, where N is the number of
%   data points and P is the number of dimensions (variables). The
%   partition minimizes the sum of point-to-cluster-centroid Euclidean
%   distances of all clusters. The returned N x 1 vector IDX contains the
%   cluster indices of each point.
%
%   IDX = k_means(X, C) works with the initial centroids, C, (K x P).
%
%   [IDX, C] = k_means(X, K) also returns the K cluster centroid locations
%   in the K x P matrix, C.
%
% See also kmeans

% Version 2.0, by Yi Cao at Cranfield University on 27 March 2008.

% Example 1: small data set
%{
N=200;
X = [randn(N,2)+ones(N,2); randn(N,2)-ones(N,2)];
[cidx, ctrs] = k_means(X, 2);
plot(X(cidx==1,1),X(cidx==1,2),'r.',X(cidx==2,1),X(cidx==2,2),'b.', ctrs(:,1),ctrs(:,2),'kx');
%}

% Example 2: large data set
%{
N=20000;
X = [randn(N,2)+ones(N,2); randn(N,2)-ones(N,2)];
tic
[cidx, ctrs] = k_means(X, 2);
toc
plot(X(cidx==1,1),X(cidx==1,2),'r.',X(cidx==2,1),X(cidx==2,2),'b.', ctrs(:,1),ctrs(:,2),'kx');
%}

% Example 3: large data set with 5 centroids 
%{
N=20000;
X = [randn(N,2)+ones(N,2); randn(N,2)-ones(N,2)];
tic
[cidx, ctrs] = k_means(X, 5);
toc
plot(X(cidx==1,1),X(cidx==1,2),'.',...
X(cidx==2,1),X(cidx==2,2),'.',...
X(cidx==3,1),X(cidx==3,2),'.',...
X(cidx==4,1),X(cidx==4,2),'.',...
X(cidx==5,1),X(cidx==5,2),'.',...
ctrs(:,1),ctrs(:,2),'+','linewidth',2)
%}

% Example 4: Comparison with kmeans in Statistics Toolbox
%{
N=20000;
X = [randn(N,2)+ones(N,2); randn(N,2)-ones(N,2)];
rand('state',0);
tic
cidx = k_means(X, 20);
toc
% Compare with kmeans in Statistis Toolbox
rand('state',0);
tic,
cidx1 = kmeans(X, 20, 'Option', statset('MaxIter',200));
toc
%}

% Check input and output
narginchk(2,2);

[n,m]=size(X);

% Check if second input is centroids
if ~isscalar(k)
    c=k;
    k=size(c,1);
else
    c=X(ceil(rand(k,1)*n),:);
end

% allocating variables
g0=ones(n,1);
gIdx=zeros(n,1);
D=zeros(n,k);

% Main loop converge if previous partition is the same as current
while any(g0~=gIdx)
%     disp(sum(g0~=gIdx))
    g0=gIdx;
    % Loop for each centroid
    for t=1:k
        d=zeros(n,1);
        % Loop for each dimension
        for s=1:m
            d=d+(X(:,s)-c(t,s)).^2;
        end
        D(:,t)=d;
    end
    % Partition data to closest centroids
    [~,gIdx]=min(D,[],2);
    % Update centroids using means of partitions
    for t=1:k
        c(t,:)=mean(X(gIdx==t,:));
    end
%     for t=1:m
%         c(:,t)=accumarray(gIdx,X(:,t),[],@mean);
%     end
end

end


function ix = otsu(x,n)

if mod(n,2)~=0
    error('number of classes is odd');
end

I = x >= graythresh(x);
ix = I*n/2;% + 1;
% do until all classes are 
if n > 2
ix(I) = otsu(x(I),n/2) + ix(I); 
ix(~I) = otsu(x(~I),n/2) + ix(~I);
else
    ix = ix+1;
    
end

end




