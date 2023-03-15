% Load DEM
DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
% Calculate flow directions and accumulations. Will be needed later.
FD  = FLOWobj(DEM);
A   = flowacc(FD);

% Modify the DEM so that it represents a block with a narrow, 1000 m lower
% rim around it.
DEM = DEM*0 + 1000;
DEM.Z(20:end-19,20:end-19) = 2000;

% Plot the DEM.
imageschs(DEM)
niceticks
pause(1);

% Let's start with a maximum gradient (excesstopography with Universal
% Slope)
maxs = 0.2;
DEM2 = excesstopo(DEM,maxs);
imageschs(DEM2);
pause(1);

% The same, but with a larger slope
maxs = 0.5;
DEM2 = excesstopo(DEM,maxs);
imageschs(DEM2);
pause(1);

% Now let's vary maximum slope spatially
[X,Y] = getcoordinates(DEM,'mat');
maxs = sin(X/2000) + 1;
DEM2 = excesstopo(DEM,maxs);
imageschs(DEM2)
pause(1);

% Same as before, but with generally higher slopes
maxs = maxs + 0.1;
DEM2 = excesstopo(DEM,maxs);
imageschs(DEM2)
pause(1);

% Yet another variable slope
maxs = (sin(X/2000) + cos(Y/2000) + 2)*0.2;
DEM2 = excesstopo(DEM,maxs);
imageschs(DEM2)
pause(1);

% A slope gradient
maxs = (X-min(X(:)))*0.000003;
DEM2 = excesstopo(DEM,maxs);
imageschs(DEM2)
pause(1);

% Random maximum slopes
maxs = rand(DEM.size);
DEM2 = excesstopo(DEM,maxs);
imageschs(DEM2)
pause(1);

% Random autocorrelated loggaussian slopes
maxs = exp(imgaussfilt(randn(DEM.size),15));
maxs = (maxs-min(maxs(:)))/(max(maxs(:))-min(maxs(:)))*0.1;
DEM2 = excesstopo(DEM,maxs);
imageschs(DEM2)
pause(1);

% Maximum slopes as a function of drainage area computed from the actual
% Big Tujunga DEM
maxs  = A^(-0.5);
DEM2 = excesstopo(DEM,maxs);
imageschs(DEM2)
pause(1);

% Maximum slopes as a function of drainage area computed from the actual
% Big Tujunga DEM 
maxs  = 2*A^(-0.5) + 0.01;
DEM2 = excesstopo(DEM,maxs);
imageschs(DEM2)
pause(1);

% Maximum slopes as a function of drainage area computed from the actual
% Big Tujunga DEM + noise
maxs  = A^(-0.5) + rand(DEM.size)*0.1;
DEM2 = excesstopo(DEM,maxs);
imageschs(DEM2)
pause(1);

% Maximum slopes as a function of drainage area computed from the actual
% Big Tujunga DEM + noise
maxs  = A^(-0.5) + rand(DEM.size)*0.01;
maxs  = min(maxs,0.2);
DEM2 = excesstopo(DEM,maxs);
imageschs(DEM2)
pause(1);