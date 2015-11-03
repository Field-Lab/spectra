function splitNeurons( rootFolder, datasets, timeTags )
    % SPLITNEURONS Splits the global neurons.mat file into each dataset's own
    % during concatenated analysis
    %
    % We get the merged neurons file in the root folder
    % then filter by time intervals and save the appropriate
    % neurons.mat in each dataset folder.
    
    nDatasets = numel(datasets);
    load([rootFolder,filesep,'concat.neurons.mat']);
    
    globalSpikeTimes = neuronSpikeTimes; % copy for names overlaps
    
    % Process time intervals
    nOffset = zeros(1,nDatasets);
    for d = 1:(nDatasets-1)
        load([datasets{d},'.spikes.mat'],'nSamples');
        nOffset(d+1) = nSamples;
    end
    nOffset = cumsum(nOffset);
    nOffset(1) = -str2double(timeTags{1}(2:find(timeTags{1} == '-',1)-1)) * 20000; % HARDCODED SAMPLING RATE
    startTimes = [0,nOffset(2:end)];
    stopTimes = [nOffset(2:end),Inf];
    
    for d = 1:nDatasets
        neuronSpikesTimes =...
            cellfun(@(x) x(and(x >= startTimes(d),x < stopTimes(d))) - nOffset(d),...
            globalSpikeTimes,'uni',false);
        save([datasets{d},'.neurons.mat'],'neuronEls','neuronClusters','neuronSpikeTimes');
    end
end

