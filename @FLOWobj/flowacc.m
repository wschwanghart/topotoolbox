function OUT = flowacc(FD,W0,RR)

%FLOWACC flow accumulation (upslope area, contributing area)
%
% Syntax
%
%     A = flowacc(FD)
%     A = flowacc(FD,W0)
%     A = flowacc(FD,W0,RR)
%
% Description
%
%     flowacc calculates the number of upstream cells based on a flow
%     direction object (FLOWobj). To obtain upslope or contributing area in
%     metric units, multiply A by the square of the cellsize, e.g.
%     A = flowacc(FD).*(FD.cellsize^2).
%
%     The second input argument can be used to define spatially variable
%     weights into the flow accumulation e.g., to simulate spatially
%     variable precipitation patterns. By default, W0 is a grid of ones.
%
%     The third input argument is the runoff ratio. By default, the runoff
%     ratio equals one everywhere. To simulate infiltration or channel
%     transmission losses, values between 0 and 1 indicate the proportion
%     of flow transferred from a cell to its downstream neighbor. 
%
% Input arguments
%
%     FD    Flow direction object (class: FLOWobj)
%     W0    weight grid (class: GRIDobj) 
%     RR    runoff ratio grid (class: GRIDobj)
%
% Output arguments
%
%     A     flow accumulation grid (class: GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     A = flowacc(FD);
%     imageschs(DEM,dilate(sqrt(A),ones(5)),'colormap','flowcolor')
%     
% 
% See also: FLOWobj, GRIDobj, FLOWobj/drainagebasins
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. March, 2016
  

% 4/3/2016: the function now makes copies of FD.ix and FD.ixc (see 
% https://topotoolbox.wordpress.com/2015/10/28/good-and-possibly-bad-news-about-the-latest-matlab-r2015b-release/comment-page-1/#comment-127


% check input arguments
narginchk(1,3)
if nargin >= 2 
    if ~isempty(W0)
        validatealignment(FD,W0);
    end
    if nargin == 3 && ~isempty(RR)
        if isa(RR,'GRIDobj')
            validatealignment(FD,RR);
        end
    end
end

if ~(exist(['flowacc_mex.' mexext],'file') == 3 && nargin<3 && strcmp(FD.type,'single'))
    if nargin == 1 || (nargin > 1 && isempty(W0))
        A = ones(FD.size);
    else        
        if isa(W0,'GRIDobj')
            A = double(W0.Z);
        else            
            A = double(W0);
        end        
        
    end
    
    if nargin == 3
        if isa(RR,'GRIDobj')
            RR = RR.Z;
        end   
    end
    
    % copies of ix and ixc to increase speed with 2015b
    ix = FD.ix;
    ixc = FD.ixc;
    
    switch FD.type
        case 'single'
            if nargin < 3
                
                for r = 1:numel(ix)
                    A(ixc(r)) = A(ix(r))+A(ixc(r));
                end
            else
                if isscalar(RR)
                    % if RR is a scalar, RR is assumed to be the
                    % coefficient of a homogenous differential equation
                    % dA/dx = RR*A
                    %
                    % A1/A2 = A0exp(RR*x1)/A0exp(RR*x2)
                    %       = exp(RR*x1)/exp(RR*x2)
                    %       = exp(RR*(x1-x2))
                    %       = exp(RR*dx)
                    
                    dx = getdistance(ix,ixc,FD.size,FD.cellsize);
                    RR = exp(RR*dx);
                    clear dx
                    
                    for r = 1:numel(ix)
						A(ixc(r)) = A(ix(r))*RR(r)+A(ixc(r));
                    end
                else
                    for r = 1:numel(ix)
                        A(ixc(r)) = A(ix(r))*RR(ix(r)) + A(ixc(r));
                    end
                end
                
                
            end
            
        case {'multi','Dinf'}
            fraction = FD.fraction;
            if nargin < 3
                for r = 1:numel(ix)
                    A(ixc(r)) = A(ix(r))*fraction(r) + A(ixc(r));
                end
            else
                for r = 1:numel(ix)
                    A(ixc(r)) = A(ix(r))*fraction(r)*RR(ix(r)) + A(ixc(r));
                end
            end
    end
else
    if nargin == 2
        if isa(W0,'GRIDobj')
            A  = flowacc_mex(FD.ix,FD.ixc,double(W0.Z));
        else            
            A  = flowacc_mex(FD.ix,FD.ixc,double(W0));
        end
    elseif nargin == 1
        A  = flowacc_mex(FD.ix,FD.ixc,FD.size);
    end
end


%% Prepare Output
% empty GRIDobj
OUT = copy2GRIDobj(FD);
% write output to GRIDobj
OUT.Z = A;
OUT.zunit = 'nr of cells';
OUT.name  = 'flow accumulation';



end