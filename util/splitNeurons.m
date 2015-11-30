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
    globalNeuronEls = neuronEls;
        globalNeuronClusters = neuronClusters;

    % Process time intervals
    nOffset = zeros(1,nDatasets);
    for d = 1:(nDatasets-1)
        load([datasets{d},'.spikes.mat'],'nSamples');
        nOffset(d+1) = nSamples;
    end
    nOffset = cumsum(nOffset);
    
        
        if any(timeTags{1} == '-') && timeTags{1}(2) ~= '-'
        nOffset(1) = -str2double(timeTags{1}(2:find(timeTags{1} == '-',1)-1    )) * 20000; % HARDCODED SAMPLING RATE
    else
        nOffset(1) = 0;
    end
        
        startTimes = [0,nOffset(2:end)];
    stopTimes = [nOffset(2:end),Inf];

    for d = 1:nDatasets
        neuronSpikeTimes =...
            cellfun(@(x) x(and(x >= startTimes(d),x < stopTimes(d))) - nOffset(d),...
            globalSpikeTimes,'uni',false);
                empties = cellfun(@isempty, neuronSpikeTimes);
                neuronSpikeTimes(empties) = [];
                neuronEls = globalNeuronEls;
                neuronEls(empties) = [];
                neuronClusters = globalNeuronClusters;
                neuronClusters(empties) = [];
        save([datasets{d},'.neurons.mat'],'neuronEls','neuronClusters','neuronSpikeTimes');
    end
end

