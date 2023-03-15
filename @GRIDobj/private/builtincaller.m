function OUT = builtincaller(funname,varargin)

% function to call built-in functions 
%
% Syntax
%
%     OUT = builtincaller(funname,G1,G2)
%
% Example
%
%     OUT = builtincaller('plus',G1,G2)
%
%


isGRIDobj = cellfun(@(x) isa(x,'GRIDobj'),varargin);
ref       = find(isGRIDobj,1,'first');
OUT       = varargin{ref};

if all(isGRIDobj)
    validatealignment(varargin{1},varargin{2});
    OUT.Z = builtin(funname,varargin{1}.Z,varargin{2}.Z);
    
else
    if ~isscalar(varargin{~isGRIDobj})
        validatealignment(varargin{isGRIDobj},varargin{~isGRIDobj});
    end
    
    if ref == 1
        OUT.Z = builtin(funname,varargin{isGRIDobj}.Z,varargin{~isGRIDobj});
    else
        OUT.Z = builtin(funname,varargin{~isGRIDobj},varargin{isGRIDobj}.Z);
    end
end

OUT.name = funname;
OUT.zunit = '';
