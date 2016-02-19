function noiseSpikes = generateNoiseEvents( realSpikes )
    %GENERATENOISEEVENTS Generates random spike times on all electrodes
    %   and makes sure they are not too close to real spikes
    %
    % Author -- Vincent Deo -- Stanford University -- August 27, 2015
    
    % Load config
    global GLOBAL_CONFIG
    covConfig = GLOBAL_CONFIG.getCovConfig();
    
    electrodes = unique(realSpikes(:,2),'sorted')';
    
    noiseSpikes = cell(max(electrodes),1);
    
    % spikes are supposed to be in order - no sorting required to get time edges
    minTime = realSpikes(1,1);
    maxTime = realSpikes(end,1);
    
    % Prorated number of noise events to generate vs a 1/2 hour run
    nNoise = round(covConfig.noiseEvents * (maxTime - minTime + 1) / (1800 * 20000));
    
    for el = electrodes
        if el == 1
            continue
        end
        noiseSpikes{el} = [randsample(maxTime-minTime, nNoise) + minTime,...
            el .* ones(nNoise,1)]; % Random generation
        
        % Round off to covConfig.noiseSpacing samples
        spikesElReduced = round(realSpikes(realSpikes(:,2) == el,1) / covConfig.noiseSpacing);
        
        % Remove duplicates
        [nRed,noiseInd,~] = unique(round(noiseSpikes{el}(:,1) / covConfig.noiseSpacing));
        spRed = unique(spikesElReduced);
        
        % Remove noise spikes overlapping real spikes
        [~,~,in] = intersect(spRed,nRed);
        noiseSpikes{el}(noiseInd(in),:) = [];
    end % el
    
    noiseSpikes = vertcat(noiseSpikes{:}); % Merge electrodes
    if numel(noiseSpikes) > 0
        noiseSpikes = sortrows(noiseSpikes,1); % Required sorting for Covariance calculation
    else
        noiseSpikes = zeros(0,2);
    end
end

