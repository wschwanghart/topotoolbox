function dem = excesstopoUniSlope(dem,cs,maxs,varargin)

if nargin == 3
    maxiter = 1000;
else
    maxiter = varargin{1};
end

% neighbor offsets (clockwise from upper left neighbor)
% row offset
neighr  = [-1 -1 -1 0 1 1 1 0];
% col offset
neighc  = [-1 0 1 1 1 0 -1 -1];
csd     = sqrt(2*cs^2);
% Diagonal maximum rise
mxrd    = csd*maxs;
% Cardinal maximum rise
mxrc    = cs*maxs;
maxrise = repmat([mxrd mxrc],1,4);

dem_ini = dem;

for iter = 1:maxiter

    % sort by elevation
    [~,ix] = sort(dem(:),'ascend');

    [rowmax,colmax] = size(dem);

    rqueue = [];
    cqueue = [];
    % rqueue = zeros(numel(dem),1);
    % cqueue = zeros(numel(dem),1);

    INQUEUE = false(size(dem));

    % Loop starts here
    r  = 1;
    r2 = 1;
    % fillqueue_counter = 1;

    while r < numel(ix)

        if r2 > numel(rqueue)
            % If the second counter has worked through the second queue
            if ~isempty(rqueue)
                INQUEUE(fastsub2ind(rqueue,cqueue,rowmax)) = false;
            end

            % reset queue
            rqueue = [];
            cqueue = [];
            r2 = 1;

            % current pixel
            ixc = ix(r);

            % increase counter
            r = r + 1;

            % Get subscripts of ixc
            [cr,cc] = fastind2sub(ixc,rowmax);

        else

            cr = rqueue(r2);
            cc = cqueue(r2);

            % ISMIN(cr,cc) = true;

            r2 = r2+1;

        end

        % Get neighbor subscripts
        nr = cr + neighr; % neighbor row
        nc = cc + neighc; % neighbor column

        for niter = 1:8

            % First check whether neighbors are outside the grid
            % Skip neighbor if true
            if nr(niter) < 1 || nr(niter) > rowmax || ...
                    nc(niter) < 1 || nc(niter) > colmax
                continue
            end

            % Elevation of central pixel
            zc = dem(cr,cc);
            % Elevation of neighbor pixel
            zn = dem(nr(niter),nc(niter));
            % Elevation of central pixel + maximum rise
            znmax = zc+maxrise(niter);

            % Needs the neighbor pixel to be lowered?
            if znmax >= zn
                % no lowering required
                continue
            end

            % If the neighbor pixel is higher
            dem(nr(niter),nc(niter)) = znmax;

            % Is the neighbor pixel now lower than the next pixel in the main
            % queue? If yes, add the pixel to the interior queue
            if dem(nr(niter),nc(niter)) <= dem(ix(r))

                if ~INQUEUE(nr(niter),nc(niter))
                    INQUEUE(nr(niter),nc(niter)) = true;
                    % add to queue
                    rqueue = [rqueue; nr(niter)];
                    cqueue = [cqueue; nc(niter)];
                end

            end

        end

    end

    if isequal(dem,dem_ini)
        break
    else
        dem_ini = dem;
    end
end


end

function [r,c] = fastind2sub(ixc,rowmax)
r = rem(ixc,rowmax);
if r == 0
    r = rowmax;
end
c = ceil(ixc/rowmax);
end

function ix = fastsub2ind(r,c,rowmax)
ix = r + (c - 1).*rowmax(1);

end