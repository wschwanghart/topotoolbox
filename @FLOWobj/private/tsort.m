function varargout = tsort(M)

% Topological sorting of a directed graph
%
% Syntax
%
%     p = tsort(M)
%     [ix,ixc] = tsort(M)
%     [ix,ixc,m] = tsort(M)
%
% Description
%
%     Topological sorting refers to a linear ordering of the vertices of
%     the directed graph of the flow direction matrix M such that each
%     upstream node becomes for downstream nodes in the ordering.
%
% Input
%
%     M       flow direction matrix
%
% Output
%
%     p       reordering
%     ix      cell index
%     ixc     cell child index
%     m       proportion 
%
%
% Example
%
%     load exampleDEM
%     M = flowdir(X,Y,dem);
%     p = tsort(M);
%     subplot(1,2,1);
%     spy(M);
%     title('M')
%     subplot(1,2,2);
%     spy(M(p,p))
%     title('M(p,p)')
%
%     pp = zeros(size(dem));
%     pp(p) = 1:numel(dem);
%
% See also: FLOWDIR
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 1. June, 2012

nrc = size(M,1);
% topological sorting of a DAG

% find vertices that have no parents
nrparent   = full(sum(spones(M),1))';
ixnoparent = find(nrparent==0);
ixpwrite   = numel(ixnoparent);

% preallocate permutation vector
p = zeros(nrc,1);

% write first noparent cells into p
p(1:ixpwrite) = ixnoparent;

% create vector to index into vertices
if nargout <= 2;
    [ixchild,ix] = find(M');
else 
    [ixchild,ix,m] = find(M');
end

nredges      = numel(ix);
ixix         = zeros(nrc,1);

% basically the same as asmember, but a little faster since the input
% vectors are already sorted. In addition the first location of occurrence 
% is returned 
t = 0;
for r = 1:numel(ix);
    if ix(r)>t;
        ixix(ix(r)) = r;
        t = ix(r);
    end
end

% Visit all childs of parentless nodes and add them to p as soon as they
% are parentless, too.
% index into p
ixpread   = 1;
ixpwrite  = ixpwrite+1;
keepgoing = true;
while keepgoing;
    
    % identify childs current parentless index
    if ixix(p(ixpread))==0;
        % node has no childs, do nothing
    else
        % ixt -> index in ix and ixchild
        ixt = ixix(p(ixpread));
        keepgoing2 = true;
        while keepgoing2
            nrparent(ixchild(ixt)) = nrparent(ixchild(ixt)) - 1;
            if nrparent(ixchild(ixt)) <= 0;
                p(ixpwrite) = ixchild(ixt);
                ixpwrite = ixpwrite + 1;
            end
            
            ixt = ixt+1;
            if ixt>nredges
                keepgoing2 = false;
            else
                if ix(ixt) ~= ix(ixix(p(ixpread)))
                    keepgoing2 = false;
                end
            end
        end
    end
    ixpread = ixpread + 1;
    
    if ixpread>nrc
        keepgoing = false;
    end
end


% prepare output
if nargout == 1;
    varargout{1} = p;
else
    clear ix ixchild
    M = M(p,:);
    if nargout == 2
        [ixchild,ixp] = find(M');
    else
        [ixchild,ixp,m] = find(M');
    end
    
    ixp = p(ixp);
    
    varargout{1} = ixp;
    varargout{2} = ixchild;
    
    if nargout == 3;
        varargout{3} = m;
    end
end
    
        
    
    
    
    
    
        

