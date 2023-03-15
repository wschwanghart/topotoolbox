function IX = coord2ind(FD,x,y)

[X,Y] = refmat2XY(FD.refmat,FD.size);
IX = coord2ind(X,Y,x,y);