function ixc = findsteepestneighbor(dem,D,cs,ix)


[nrow,ncol] = size(dem);
siz = [nrow,ncol];
[row,col]   = ind2sub(size(dem),ix);

coln = [0, 1, 1, 1, 0, -1, -1, -1];
rown = [-1, -1, 0, 1, 1, 1, 0, -1];
csd  = hypot(cs,cs);
dist = [cs csd cs csd cs csd cs csd]; 

ixc  = zeros(size(ix),class(ix));

for ixix = 1:numel(ix);
    if isinf(D(ix(ixix)));
        inflat = false;
        g = 0;
    else
        inflat = true;
        g = 0;
    end
    
        for neighbor = 1:8;
            c = col(ixix) + coln(neighbor);
            if c < 1 || c > ncol;
                continue
            end
            r = row(ixix) + rown(neighbor);
            if r < 1 || r > nrow;
                continue
            end
            
            ixn = (c-1)*nrow + r;
%             ixn = sub2ind(siz,r,c);
            % if your index is in nonflat areas
            if ~inflat;
                gg = (dem(ix(ixix))-dem(ixn))/dist(neighbor);
                if gg > g;
                    ixc(ixix) = ixn;
                    g = gg;
                end
            else
                if D(ix(ixix)) == 1;
                    if dem(ixn) == dem(ix(ixix)) && isinf(D(ixn));
                        ixc(ixix) = ixn;
                        break
                    end
                else
                
                gg = (D(ix(ixix))-D(ixn))/dist(neighbor);%D(ixn);
                if gg > g && ~isinf(gg);
                    ixc(ixix) = ixn;
                    g = gg;
                end
                end
            end
            
        end
end

            
            
        
        

