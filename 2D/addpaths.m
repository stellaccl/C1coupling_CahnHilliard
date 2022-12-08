% add all needed function paths
try
    functionname='addpaths.m';
    functiondir=which(functionname);
    functiondir=functiondir(1:end-length(functionname));
    addpath([functiondir '/PHTutils'])
    addpath([functiondir, '/inst'])
    
catch end