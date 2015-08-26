function noiseSpikes = generateNoiseEvents( realSpikes )
    %GENERATENOISEEVENTS Generates random spike times on all electrodes
    %   and makes sure they are not too close to real spikes
    
    config = mVisionConfig();
    covConfig = config.getCovConfig();
    
    electrodes = unique(realSpikes(:,2),'sorted')';
    
    noiseSpikes = cell(max(electrodes),1);
    
    % spikes are supposed to be in order
    minTime = realSpikes(1,1);
    maxTime = realSpikes(end,1);
    
    % Prorated vs a 1/2 hour run
    nNoise = round(covConfig.noiseEvents * (maxTime - minTime + 1) / 36000000);
    
    for el = electrodes
        if el == 1
            continue
        end
        noiseSpikes{el} = [randsample(maxTime-minTime, nNoise) + minTime,...
            el .* ones(nNoise,1)];
        spikesElReduced = round(realSpikes(realSpikes(:,2) == el,1) / 20); % 20 so that any noise spike within 10 samples of a real spike will be discarded
        
        [nRed,noiseInd,~] = unique(round(noiseSpikes{el}(:,1) / 20));
        spRed = unique(spikesElReduced);
        [~,~,in] = intersect(spRed,nRed);
        noiseSpikes{el}(noiseInd(in),:) = [];
    end
    
    noiseSpikes = vertcat(noiseSpikes{:});
    noiseSpikes = sortrows(noiseSpikes,1);
    
end

