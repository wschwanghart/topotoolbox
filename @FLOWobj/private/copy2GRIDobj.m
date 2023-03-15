function G = copy2GRIDobj(FD)

% copy common properties from FLOWobj to GRIDobj

% empty GRIDobj
G = GRIDobj([]);
% find common properties of F and G and from F to G
pg = properties(G);
pf = properties(FD);
p  = intersect(pg,pf);
for r = 1:numel(p)
    G.(p{r}) = FD.(p{r});
end