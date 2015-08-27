function projSpikes = loadPrj( prjFilePath, el )
    %LOADPRJ Loads the PC projections for a given electrode in a given file
    %   Purpose of this wrapper is to maintain transparency within parfor loops,
    %   If one is used to perform parallelization over electrodes in PCClustering.m
    %
    % Author -- Vincent Deo -- Stanford University -- August 21, 2015
    
    eval(['load(''', prjFilePath,sprintf(''',''projSpikes%u''',el),');']);
    eval(sprintf('projSpikes = projSpikes%u;',el));
end

