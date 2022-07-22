function tthelp(str,page)
%TTHELP Search on the TopoToolbox blog
%
% Syntax
%
%     tthelp keyword
%     tthelp('keyword')
%
% Description
%
%     TTHELP lets you search the TopoToolbox blog and returns links to the
%     blog to the command line. 
%
% Example
%
%     tthelp ksn
%
% Note that this function may not work properly if wordpress does changes
% to its html-code.
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 22. July, 2022

if nargin == 1
    page = 1;
end

% First page URL
if page == 1
    url = 'https://topotoolbox.wordpress.com/?s=';
else
    pagestr = num2str(page);
    url = ['https://topotoolbox.wordpress.com/page/' pagestr '/?s='];
end

t = webread([url str]);

[tokens,matches] = regexp(t,'<h3>(.*?)</h3>','tokens','match');
for r = 1:numel(tokens)
    % replace &nbsp; with blank
    tokens{r} = replace(tokens{r},'&nbsp;',' ');
    disp(tokens{r}{1})
end

% now check whether there's another page with results
findnextstr = '<li><a class="next page-numbers" (.*?)>Next &rarr;</a></li>';

if ~isempty(regexp(t,findnextstr,'once'))
    tthelp(str,page+1)
end



