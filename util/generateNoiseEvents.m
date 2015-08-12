function noiseSpikes = generateNoiseEvents( realSpikes )
    %GENERATENOISEEVENTS Generates random spike times on all electrodes
    %   and makes sure they are not too close to real spikes
    
    config = mVisionConfig();
    covConfig = config.getCovConfig();
    
    electrodes = unique(realSpikes(:,2),'sorted')';
    
    noiseSpikes = cell(max(electrodes),1);
    
    for el = electrodes
        if el == 1
            continue
        end
        noiseSpikes{el} = [randsample(max(realSpikes(:,1))-min(realSpikes(:,1)),...
            covConfig.noiseEvents) + min(realSpikes(:,1)),...
            el .* ones(covConfig.noiseEvents,1)];
    end
    
    noiseSpikes = vertcat(noiseSpikes{:});
    noiseSpikes = sortrows(noiseSpikes,1);
    
end

