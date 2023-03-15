function Lgrid = shufflelabel(Lgrid,r)

%SHUFFLELABEL shufflelabel randomly relabels a label matrix
%
% Syntax
%
%     L = shufflelabel(L)
%     L = shufflelabel(L,reset)
%
% Description
%
%     shufflelabel randomly changes the order of labels in the  
%     label matrix L. Zeros and nans are ignored. L can be every
%     numerical data type, char or a cell array of strings.
%
%     When called with two input arguments, reset is either true or 
%     false. If true, label values in L are reset to range from
%     one to numel(unique(L(:)). Note that this only works for 
%     numeric arrays.
%
% Input arguments
%
%     L      grid with ordinal or categorical data (class: GRIDobj)
%     reset  false (default) or true
%
% Output arguments
%
%     Ls     grid with values randomly shuffled
%
% Example
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     L   = reclassify(DEM,'equalquantiles',10);
%     L   = shufflelabel(L);
%     imagesc(L)
%
%
% See also: RANDPERM, BWLABEL, LABELMATRIX
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. August, 2017



% essentially, error checking is not necessary
narginchk(1,2)

% get values in GRIDobj
L = Lgrid.Z;

% is L a numeric array? 
inum = isnumeric(L);

% reset labeling?
if nargin == 1;
    % by default no.
    r = false;
else
    % applies only if L is a numeric array.
    validateattributes(r,{'numeric','logical'},{'scalar'})
    r = r && inum;
end

% size of L
siz = size(L);

% force column vector
L   = L(:);

% exclude zeros and nans from shuffling (only when L is 
% a numeric array)
if inum
    I   = ~(L==0 | isnan(L));
else
    I   = ':';
end

% find unique elements in label vector
[uniqueL,~,ix] = unique(L(I));

% shuffle labels
if r
    uniqueLS = randperm(numel(uniqueL));
else
    uniqueLS = uniqueL(randperm(numel(uniqueL)));
end

% and map labels back into L
L(I) = uniqueLS(ix);

% finally reshape back to original size
L = reshape(L,siz);

% write to GRIDobj
Lgrid.Z = L;

