function [covMatrix,averages,totSpikes] = buildCovariances(spikesTotal, dataPath, timeCommand)
    % Build the covariance matrix for spikes around a given electrode
    
    %% Argument validation
    if ~(exist(dataPath,'file') == 2 || exist(dataPath,'file') == 7)
        throw(MException('','CovarianceCalculation: data folder|file does not exist'));
    end
    
    %% Load covariance configuration
    config = mVisionConfig();
    covConfig = config.getCovConfig();
    
    rawDataSource = [dataPath,timeCommand];
    
    nLPoints = covConfig.nLPoints;
    nRPoints = covConfig.nRPoints;
    nPoints = nLPoints + nRPoints + 1;
    
    
    %% Creating data source
    dataSource = DataFileUpsampler(rawDataSource, covConfig.meanTimeConstant, nLPoints, nRPoints);
        
    nElectrodes = dataSource.nElectrodes;
    disconnected = dataSource.disconnected;
    
    aligner = spikeAligner(dataSource);
    
    %% Setting up neighbor map
    % Subfunction encapsulates java use
    [adjacent,~] = catchAdjWJava( dataSource, covConfig.electrodeUsage);
    
    %% Initializing covariance matrices
    covMatrix = cell(nElectrodes,1);
    averages = cell(nElectrodes,1);
    totSpikes = zeros(nElectrodes,1);
        
    for el = 2:nElectrodes
        if disconnected(el)
            continue
        end
        
        covMatrix{el} = zeros((nPoints - 2) * numel(adjacent{el}));
        averages{el} = zeros(1,(nPoints - 2) * numel(adjacent{el}));
    end
    
    while ~dataSource.isFinished % stopSample should be the first sample not loaded
        
        [bufferStart,bufferEnd] = dataSource.loadNextBuffer();
        
        dataSource.createInterpolant();
        
        %% Load Spikes
        spikesTemp = spikesTotal(and(spikesTotal(:,1) >= bufferStart + nLPoints, spikesTotal(:,1) < bufferEnd - nRPoints),:);
        % If no spikes at all are loaded, skip buffer
        if numel(spikesTemp) == 0
            continue
        end
        
        spikes = cell(nElectrodes,1);
        
        for el = 1:nElectrodes
            spikes{el} = spikesTemp(spikesTemp(:,2) == el,1);
        end
        
        clear spikesTemp;
        
        %% Process by electrodes
        % Could parallel here, but actually slower due to the IO cost of sending to each worker.
        % The better part would be to change the while loop to a smarter parfor loop.
        for el = 2:nElectrodes
            
            nSpikes = numel(spikes{el});
            
            if disconnected(el) || nSpikes == 0
                continue
            end            
            
            %% Align all spikes with the aligner            
            spikesTemp = aligner.alignSpikes(el,spikes{el});
            
            totSpikes(el) = totSpikes(el) + nSpikes;
            
            averages{el} = averages{el} + sum(spikesTemp,1);
            covMatrix{el} = covMatrix{el} + spikesTemp' * spikesTemp;
            
        end % el
    end % while ~isFinished
    
    %% Normalize averages and covmatrix
    
    for el = 2:nElectrodes
        if dataSource.disconnected(el) || totSpikes(el) == 0
            continue;
        end
        
        if totSpikes(el) >= 2
            covMatrix{el} = (covMatrix{el} - averages{el}' * averages{el} / totSpikes(el))/(totSpikes(el)-1);
            averages{el} = averages{el}/totSpikes(el);
        end
    end
end