function thresholds = multilevel_otsu(h,nClasses,values)

%MULTILEVEL_OTSU Multilevel Otsu thresholding
%
% Syntax
%
%     thresholds = multilevel_otsu(H,nClasses)
%     thresholds = multilevel_otsu(H,nClasses,values)
%
% Description
%
%     multilevel_otsu computes thresholds for a multilevel classification
%     of the histogram. nClasses is the number of classes and must be an
%     integer value between 2 and 5.
%
%     By default, the function assumes that values of the histogram are
%     integers from 1 to numel(H). The vector values can contain different
%     values, but must have uniform spacing. 
%
% Input arguments
%
%     H         Vector with absolute or relative frequencies
%     nClasses  Number of classes (2-5)
%
% Output arguments
%
%     thresholds  Vector with thresholds with nClasses-1 elements.
%
% Reference
% 
%     The function has been adopted from David Leglands matlab-image-class
%     https://github.com/mattools/matlab-image-class
%
%
% 
% See also: GRIDobj/reclassify
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 27. March, 2023

%% Check input parameters
validateattributes(nClasses,"numeric",{'scalar','>=',2,'<=',5},2);
if nargin == 3
    assert(numel(h)==numel(values),...
        'H and values must have the same number of elements.')
end

%% Initialisations

% compute histogram, and convert into probability density
h = h / sum(h);

% number of gray levels
nLevels = length(h);

Puv = zeros(nLevels, nLevels);
Suv = zeros(nLevels, nLevels);
Huv = zeros(nLevels, nLevels);

% initialize diagonal terms
for i = 1:nLevels
    Puv(i, i) = h(i);
    Suv(i, i) = h(i) * i;
end

% initialize first rows
for i = 2:nLevels
    Puv(1, i) = Puv(1, i-1) + h(i);
    Suv(1, i) = Suv(1, i-1) + i*h(i);  
end

% initialize the rest of the matrices
for u = 2:nLevels
    for v = u+1:nLevels
        Puv(u, v) = Puv(1, v) - Puv(1, u-1);
        Suv(u, v) = Suv(1, v) - Suv(1, u-1);
    end
end

% now calculate Huv
for u = 1:nLevels
    for v = u+1:nLevels
        if Puv(u,v) > 0
            Huv(u,v) = Suv(u,v) * Suv(u,v) / Puv(u,v);
        end
    end
end


%% Main processing

% create array for threshold values (first one is zeros, last one is inf)
thresholds = zeros(nClasses+1, 1);
thresholds(end) = inf;

% iterate over all possible combinations
maxSig = 0.0;
switch nClasses
    case 2
        % two classes -> only need to find second threshold value
        for i = 1:nLevels-1
            Sq = Huv(1,i) + Huv(i+1, end);
            if Sq >= maxSig
                thresholds(2) = i;
                maxSig = Sq;
            end
        end
        
    case 3
        % three classes
        for i = 1:nLevels-2
            for j = i+1:nLevels-1
                Sq = Huv(1,i) + Huv(i+1, j) + Huv(j+1, end);
                if Sq >= maxSig
                    thresholds(2) = i;
                    thresholds(3) = j;
                    maxSig = Sq;
                end
            end
        end
        
    case 4
        % four classes
        for i = 1:nLevels-3
            for j = i+1:nLevels-2
                for k = j+1:nLevels-1
                    Sq = Huv(1,i) + Huv(i+1, j) + Huv(j+1, k) + Huv(k+1, end);
                    if Sq >= maxSig
                        thresholds(2) = i;
                        thresholds(3) = j;
                        thresholds(4) = k;
                        maxSig = Sq;
                    end
                end
            end
        end
        
    case 5
        % five classes
        for i = 1:nLevels-3
            for j = i+1:nLevels-2
                for k = j+1:nLevels-1
                    for m = k+1:nLevels-1
                        Sq = Huv(1,i) + Huv(i+1, j) + Huv(j+1, k) + Huv(k+1, m) + Huv(m+1, end);
                        if Sq >= maxSig
                            thresholds(2) = i;
                            thresholds(3) = j;
                            thresholds(4) = k;
                            thresholds(5) = m;
                            maxSig = Sq;
                        end
                    end
                end
            end
        end
        
    otherwise
        error('Can not manage %d number of classes', nClasses);
end

thresholds = thresholds(2:end-1);
if nargin == 3
    thresholds = values(thresholds);
end

