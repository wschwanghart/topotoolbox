function showmethods(classname,showlink)

% displays class method names and H1 lines in the command line
%
% Syntax
%
%     showmethods(classname)
%     showmethods(classname,showlink)
%
% Description
%
%     showmethods outputs a list of methods available for a specific class
%     (e.g. GRIDobj) in the command window. Only the methods are listed
%     that feature an H1 line, e.g., the first line of the help text block
%     in a function.
%
% Input arguments
%
%     classname     string with class name (e.g. GRIDobj, FLOWobj)
%     showlink      false or true (default). If true, the function displays
%                   hyperlinks to the documentation.
%                   
%
% Example
%
%     showmethods('GRIDobj')
%
% See also: methods, properties, class
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. December, 2017

narginchk(1,2)

if nargin == 1
    showlink = true;
end


m = methods(classname);

C = cellfun(@(x) numel(x),m,'uniformoutput',true);

maxmethcharacter = max(C);

for r = 1:numel(m)
    s = which([classname '/' m{r}]);
    fileID = fopen(s,'r');
    H1line = false;
    
    try
        while ~H1line && ~feof(fileID)
            tline = fgetl(fileID);
            if isempty(tline)
                continue
            end
            if strcmp(tline(1),'%')
                H1line = true;
            end
        end
        
        if feof(fileID)
            fclose(fileID);
            continue
        end
        
        methodstr = m{r};
        if numel(methodstr)<= maxmethcharacter
            addblanks = repmat(' ',1,maxmethcharacter - numel(methodstr));
        end
        
        h1str = tline(2:end);
        ix    = strfind(h1str,' ');
        h1str = h1str(ix(1)+1:end);
        if ~showlink                
            disp([ upper(methodstr) addblanks ' : ' h1str]);
        else
            disp(['<a href="matlab: doc ' classname '/' methodstr '">' upper(methodstr) '</a>' addblanks ' : ' h1str]);
        end
        fclose(fileID);
        
    catch
        disp(m{r})
    end
    
    
end
            
           


