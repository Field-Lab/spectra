function mergeCovAndSpikes( rootFolder, datasets, timeTags)
    %MERGECOVANDSPIKES merges spike files and covariance files across datasets
    %
    % In concatenated mode, spikes and cov are computed for individual dataset,
    % and this files computes the global spikes and cov
    % by concatenation / weighted sum
    % Also, the individual cov files are overwritten with the global one,
    % so that individual prj computations end being done on identical eigenvectors.
    
    
    nDatasets = numel(datasets);

    globalSpikes = cell(1,nDatasets);
    globalTTl = cell(1,nDatasets);
    
    load([datasets{1},'.spikes.mat']);
    load([datasets{1},'.cov.mat']);
    
    globalAverages = cellfun(@(M,n) n*M,averages,num2cell(totSpikes,2),'uni',false);
    globalCovMat = cellfun(@(M,n) n*M,covMatrix,num2cell(totSpikes,2),'uni',false);
    globalTotSpikes = totSpikes;
    
    initialOffset = str2double(timeTags{1}(2:find(timeTags{1} == '-',1)-1)) * 20000; % HARDCODED SAMPLING RATE
    
    globalSpikes{1} = bsxfun(@minus,spikeSave,int32([initialOffset,0]));
    globalTTl{1} = ttlTimes;
    
    sampleOffset = nSamples;
    
    % increment, stack matrices mul by nspikes
    for d = 2:nDatasets
        load([datasets{d},'.spikes.mat']);
        load([datasets{d},'.cov.mat']);
        
        globalAverages = cellfun(@(C,M,n) C + n*M,globalAverages,averages,num2cell(totSpikes,2),'uni',false);
        globalCovMat = cellfun(@(C,M,n) C + n*M,globalCovMat,covMatrix,num2cell(totSpikes,2),'uni',false);
        globalTotSpikes = globalTotSpikes + totSpikes;
        
        globalSpikes{d} = bsxfun(@plus,spikeSave,int32([sampleOffset 0]));
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
    
    save([rootFolder,filesep,'concat.spikes.mat'],'spikeSave','ttlTimes','nSamples','-v7.3');
    save([rootFolder,filesep,'concat.cov.mat'],'averages','covMatrix','totSpikes');
		% Put a copy of the cov file in each of the subfolders
    for d = 1:nDatasets
        load([datasets{d},'.cov.mat'],'totSpikes');
        save([datasets{d},'.cov.mat'],'averages','covMatrix','totSpikes');
    end

end
