function compilemexfiles

%COMPILEMEXFILES compile mex-functions that come with TopoToolbox 2
%
% Syntax 
%
%     compilemexfiles
%
% Description
%
%     TopoToolbox 2 comes with a few mex-functions that have been
%     written in C and that must be compiled prior to usage. While all
%     functions have been written in plain m-code, the usage of
%     mex-functions enhances the speed at which some functions are
%     evaluated. compilemexfiles compiles all these files to run on your
%     system. Prior to running compilemexfiles run
%
%     mex -setup
%
%     and locate a compiler on your system.
%
% Available mex-functions
%
%     All mex-functions are located in private directories
%
%     @FLOWobj/private
%     dependencemap_mex.c
%     drainagebasins_mex.c
%     flowacc_mex.c
%     steepestneighbor_mex.c
%     tsort_mex.c
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 30. January, 2013



location = which('compilemexfiles');
[pathstr,~,~] = fileparts(location);

newFolder = [pathstr filesep '@FLOWobj' filesep 'private']; 
oldFolder = cd(newFolder);

files = dir('*.c');

try
    for r = 1:numel(files)
        funname = files(r).name;
        mex('-largeArrayDims',funname)
        % mex('-R2018a',funname)
    end
catch err
    cd(oldFolder);
    rethrow(err)
end
cd(oldFolder);






