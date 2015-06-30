function [covMatrix,averages,totSpikes] = buildCovariances(spikesTotal, dataPath, timeCommand)
    % Build the covariance matrix for spikes around a given electrode
    
    
    %% Imports
    import edu.ucsc.neurobiology.vision.electrodemap.*
    import edu.ucsc.neurobiology.vision.io.*
%     import java.io.*
    
    %% Argument validation
    if ~(exist(dataPath,'file') == 2 || exist(dataPath,'file') == 7)
        throw(MException('','CovarianceCalculation: data folder|file does not exist'));
    end
    
    %% Load covariance configuration
    config = mVisionConfig();
    covConfig = config.getCovConfig();
    
    rawDataSource = [dataPath,timeCommand];
    
    meanTimeConstant = covConfig.meanTimeConstant;
    
    nLPoints = covConfig.nLPoints;
    nRPoints = covConfig.nRPoints;
    nPoints = nLPoints + nRPoints + 1;
    
    electrodeUsage = covConfig.electrodeUsage;
    
    %% Creating data source
    dataSource = DataFileUpsampler(rawDataSource, meanTimeConstant, nLPoints, nRPoints);
        
    nElectrodes = dataSource.nElectrodes;
    disconnected = dataSource.disconnected;
    
    %% Setting up neighbor map
    % Still need to go through java to set neighbors
    header = dataSource.rawDataFile.getHeader();
    packedArrayID = int32(header.getArrayID());
    electrodeMap = ElectrodeMapFactory.getElectrodeMap(packedArrayID);
    
    adjacent = cell(nElectrodes,1);
    maxAdjacent = 0;
    for el = 0:(nElectrodes-1)
        adjacent{el+1} = electrodeMap.getAdjacentsTo(el, electrodeUsage);
        if numel(adjacent{el+1}) > maxAdjacent
            maxAdjacent = numel(adjacent{el+1});
        end
    end
    
    %% Data flow
    
    upSampRatio = dataSource.upSampleRatio;
    upSampStep = 1/upSampRatio;
    
    % interpolation bases
    resampleBase = nLPoints:upSampStep:(nLPoints+2);
    
    % Initializing covariance matrices
    covMatrix = cell(nElectrodes,1);
    averages = cell(nElectrodes,1);
    totSpikes = zeros(nElectrodes,1);
    
    spikePile = cell(nElectrodes,1);
    spikeTot = 100;
    spikeIndex = zeros(nElectrodes,1);
    
    for i = 2:nElectrodes
        covMatrix{i} = zeros((nPoints-2) * numel(adjacent{i}));
        averages{i} = zeros(1,(nPoints - 2) * numel(adjacent{i}));
        spikePile{i} = zeros(spikeTot,(nPoints - 2) * numel(adjacent{i}));
    end
    
    while ~dataSource.isFinished % stopSample should be the first sample not loaded
        
        [bufferStart,bufferEnd] = dataSource.loadNextBuffer();
        dataSource.upsampleBuffer();
        
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
            if disconnected(el)
                continue;
            end
            %% Process each spike
            for spikeTime = spikes{el}'
                % Load master spike
                interpSpike = dataSource.upSampData(el,round(upSampRatio*(resampleBase + double(spikeTime) - bufferStart - nLPoints - 1))+1);
                
                % Find minimum and compute associated resample points
                offset = (find(interpSpike == min(interpSpike),1)-1)/upSampRatio;
                interpPoints = (1:(nPoints-2)) + offset;
                
                % Assign realigned spikes, master + neighbors
                try
                    spikePile{el}(spikeIndex(el)+1,:) =...
                        reshape(dataSource.upSampData(adjacent{el}+1,...
                        round(upSampRatio*(interpPoints + double(spikeTime) -...
                        bufferStart - nLPoints - 1))+1)',...
                        size(spikePile{el},2),1);
                catch err
                    err
                    1;
                end
                spikeIndex(el) = spikeIndex(el)+1;
                
                if spikeIndex(el) == spikeTot
                    totSpikes(el) = totSpikes(el) + spikeTot;
                    averages{el} = averages{el} + spikeTot*mean(spikePile{el},1);
                    covMatrix{el} = covMatrix{el} + (spikeTot-1)*cov(spikePile{el});
                    spikeIndex(el) = 0;
                end
                
                %                 centeredSpike = dataSource.upSampData(adjacent{el}+1,...
                %                     round(upSampRatio*(interpPoints + double(spikeTime)...
                %                     - bufferStart - nLPoints - 1))+1)';
                %                 spikeAtWork = dataSource.rawData(adjacent{el}+1,...
                %                     (spikeTime-nLPoints-bufferStart+1):(spikeTime+nRPoints-bufferStart+1));
                %                 plot(1:21,spikeAtWork(1,:),'b+',resampleBase,interpSpike,'r');
                %                 hold on
                %                 plot(1:21,spikeAtWork,'b+');
                %                 plot(1:21,spikeAtWork,'k--');
                %                 plot(2:20,centeredSpike,'r-');
                %                 hold off
                
            end % spikeTime
        end % el
    end % while ~isFinished
    
    %% Normalize averages and covmatrix
    for el = 2:nElectrodes
        if disconnected(el)
            continue
        end
        
        if spikeIndex(el) > 0
            totSpikes(el) = totSpikes(el) + spikeIndex(el);
            averages{el} = averages{el} + spikeIndex(el)*mean(spikePile{el}(1:spikeIndex(el),:),1);
            if spikeIndex(el) > 1
                covMatrix{el} = covMatrix{el} + (spikeIndex(el)-1)*cov(spikePile{el}(1:spikeIndex(el),:));
            else
                covMatrix{el} = covMatrix{el} + cov(spikePile{el}(1:spikeIndex(el),:));
            end
        end
        
        
        if totSpikes(el) >= 2
            covMatrix{el} = (covMatrix{el} - averages{el}' * averages{el} / totSpikes(el))/(totSpikes(el)-1);
            averages{el} = averages{el}/totSpikes(el);
        end
    end
    
end