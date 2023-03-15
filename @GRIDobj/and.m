function OUT = and(varargin)

funname = 'and';
narginchk(2,2)
OUT = builtincaller(funname,varargin{:});