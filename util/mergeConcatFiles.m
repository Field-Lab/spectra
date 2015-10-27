function mergeConcatFiles( rootFolder, datasets)
    %MERGEFILES merges spike files and covariance files across datasets
    nDatasets = numel(datasets);

    globalSpikes = cell(1,nDatasets);
    globalTTl = cell(1,nDatasets);
    
    load([rootFolder,filesep,datasets{1},'.spikes.mat']);
    load([rootFolder,filesep,datasets{1},'.cov.mat']);
    
    globalAverages = cellfun(@(M,n) n*M,averages,num2cell(totSpikes,2),'uni',false);
    globalCovMat = cellfun(@(M,n) n*M,covMatrix,num2cell(totSpikes,2),'uni',false);
    globalTotSpikes = totSpikes;
    
    globalSpikes{1} = spikeSave;
    globalTTl{1} = ttlTimes;
    
    sampleOffset = nSamples;
    
    % increment, stack matrices mul by nspikes
    for d = 2:nDatasets
        load([rootFolder,filesep,datasets{d},'.spikes.mat']);
        load([rootFolder,filesep,datasets{d},'.cov.mat']);
        
        globalAverages = cellfun(@(C,M,n) C + n*M,globalAverages,averages,num2cell(totSpikes,2),'uni',false);
        globalCovMat = cellfun(@(C,M,n) C + n*M,globalCovMat,covMatrix,num2cell(totSpikes,2),'uni',false);
        globalTotSpikes = globalTotSpikes + totSpikes;
        
        globalSpikes{d} = spikeSave + sampleOffset;
        globalTTl{d} = ttlTimes + sampleOffset;
        
        sampleOffset = sampleOffset + nSamples;
    end
    
    % Normalize, concatenate
    averages = cellfun(@(M,n) M./n,globalAverages,num2cell(globalTotSpikes,2),'uni',false);
    covMatrix = cellfun(@(M,n) M./n,globalCovMat,num2cell(globalTotSpikes,2),'uni',false);
    totSpikes = globalTotSpikes;
    
    spikeSave = vertcat(globalSpikes{:});
    ttlTimes = vertcat(globalTTl{:});
    nSamples = sampleOffset;
    
    save([rootFolder,'.spikes.mat'],'spikeSave','ttlTimes','nSamples','-v7.3');
    save([rootFolder,'.cov.mat'],'averages','covMatrix','totSpikes');
end