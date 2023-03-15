function OUT = plus(varargin)

funname = 'plus';
narginchk(2,2)
OUT = builtincaller(funname,varargin{:});