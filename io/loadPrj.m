function projSpikes = loadPrj( prjFilePath, el )
    %LOADPRJ Loads the PC projections for a given electrode in a given file
    %   Purpose of this wrapper is to maintain transparency within parfor loops.
    
    eval(['load(''', prjFilePath,sprintf(''',''projSpikes%u''',el),');']);
    eval(sprintf('projSpikes = projSpikes%u;',el));
end

