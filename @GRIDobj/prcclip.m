function [lims,A] = prcclip(A,prc,symmetric)

%PRCCLIP percentile clipping
%
% Syntax
%
%     lims = prcclip(A,prc)
%     lims = prcclip(A,prc,symmetric)
%     [lims,Ac] = prcclip(A,prc)
%     
% Description
%
%     PRCCLIP determines the values lims that clip the values in A to range
%     between the nth percentile and 100-nth percentile. prc can be a
%     scalar or two element vector. If prc is a scalar, then lims are
%     determined by the prc'th and 100-prc'th percentile. If prc is a
%     vector, it contains the lower and upper percentile (e.g. [2 98]).
%
% Input arguments
%
%     A          GRIDobj
%     prc        percentile (scalar or two element vector)
%     symmetric  if true, lims will be symmetric around zero. prcclip will
%                determine the percentiles and then will adjust the limits
%                so that they range between 
%                [-max(abs(lims))  max(abs(lims))]
%
% Output arguments
%
%     lims       limits
%     Ac         percentile clipped GRIDobj
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     C   = curvature(DEM);
%     lims = prcclip(C,2,true);
%     clr = ttscm('vik');
%     imageschs(DEM,C,'colormap',clr,'caxis',lims);
%     
%
% See also: GRIDobj/minmaxnorm
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 14. June, 2018

if nargin == 1
    prc = 2;
    symmetric = false;
elseif nargin == 2
    symmetric = false;
end
    
    
qclip = prc/100;
if numel(qclip) == 1
    qclip = [qclip 1-qclip];
end

qclip = sort(qclip);

I = ~isnan(A.Z);

[n,edges] = histcounts(A.Z(I(:)),'Normalization','cdf');
lval = edges(find(n>=qclip(1),1,'first'));
uval = edges(find(n<(qclip(2)),1,'last'));
if lval == uval
    warning('TopoToolbox:imageschs','percent clip returns flat matrix');

else
    if nargout == 2
        A.Z(I) = max(A.Z(I),lval);
        A.Z(I) = min(A.Z(I),uval);
    end
    
end

lims = [lval,uval];

if symmetric
    lims = max(abs(lims));
    lims = [-lims lims];
end

% if nargout == 2
%     if flatmatrix
%         A(:,:) = lval;
%     else
%         A = max(A,lims(1));
%         A = min(A,lims(2));
%     end
% end
    
    