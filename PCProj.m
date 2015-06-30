function [projSpikes,eigenValues,eigenVectors,spikeTimes] = PCProj(dataPath, timeCommand, spikesTotal, covMatrix, averages, totSpikes)
    % Build the covariance matrix for spikes around a given electrode
    % Input HashMap parameters should be the same given than for SpikeFindingM
    
    
    %% Imports
    import edu.ucsc.neurobiology.vision.electrodemap.*
    import edu.ucsc.neurobiology.vision.io.*
%     import java.io.*
    
    %% Argument validation
    if ~(exist(dataPath,'file') == 2 || exist(dataPath,'file') == 7)
        throw(MException('','CovarianceCalculation: data folder|file does not exist'));
    end
    
    %% Load projections configuration
    config = mVisionConfig();
    projConfig = config.getProjConfig();
    
    rawDataSource = [dataPath,timeCommand];
    
    meanTimeConstant = projConfig.meanTimeConstant;
    
    nLPoints = projConfig.nLPoints;
    nRPoints = projConfig.nRPoints;
    nPoints = nLPoints + nRPoints + 1;
    
    electrodeUsage = projConfig.electrodeUsage;
    
    nDims = projConfig.nDims;
    
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
    
    % Projections storage and init
    eigenValues = cell(nElectrodes,1);
    eigenVectors = cell(nElectrodes,1);
    projSpikes = cell(nElectrodes,1);
    spikeTimes = cell(nElectrodes,1);
    currSpike = ones(nElectrodes,1);
    
    for el = 2:nElectrodes
        if disconnected(el)
            continue
        end
        
        [v,d] = eig(covMatrix{el});
        e = flipud(diag(d));
        eigenValues{el} = e(1:nDims);
        eigenVectors{el} = fliplr(v(:,(end-nDims+1):end));
        projSpikes{el} = zeros(totSpikes(el),nDims);
        spikeTimes{el} = zeros(1,totSpikes(el));
    end
    
    while ~dataSource.isFinished % stopSample should be the first sample not loaded
        
        [bufferStart,bufferEnd] = dataSource.loadNextBuffer();
        dataSource.upsampleBuffer();
        
        %% Load Spikes
        spikesTemp = spikesTotal(and(spikesTotal(:,1) >= bufferStart + nLPoints, spikesTotal(:,1) < bufferEnd - nRPoints),:);
        % If no spikes at all are loaded, skip iteration
        % Required as by Matlab cast spikes is empty 513x0 and not a cell array in that case
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
                continue
            end
            %% Store spike times
            spikeTimes{el}(currSpike(el):(currSpike(el)+numel(spikes{el})-1)) = spikes{el};
            
            %% Process each spike
            for spikeTime = spikes{el}'
                % Load master spike
                interpSpike = dataSource.upSampData(el,...
                    round(upSampRatio*(resampleBase + double(spikeTime)...
                    - bufferStart - nLPoints - 1))+1);
                
                % Find minimum and compute associated resample points
                [~,offset] = min(interpSpike);
                offset = (offset-1)/upSampRatio;
                interpPoints = (1:(nPoints-2)) + offset;
                
                % Load realigned spikes, master + neighbors
                centeredSpike = dataSource.upSampData(adjacent{el}+1,...
                    round(upSampRatio*(interpPoints +...
                    double(spikeTime) - bufferStart - nLPoints - 1))+1)';
                
                % Compute projections
                projSpikes{el}(currSpike(el),:) = (centeredSpike(:)' - averages{el}) * eigenVectors{el};
                currSpike(el) = currSpike(el)+1;
                % TODO add error cases to check number of spikes - possibly not the same as stored,
                % etc. Possible array extension to implement, etc.
                
                % TODO: Check what happens for 30 first spikes?
                
                %                 spikeAtWork = rawData(adjacent{el}+1,...
                %                     (spikeTime-nLPoints-bufferStart+1):(spikeTime+nRPoints-bufferStart+1));
                %                 plot(1:21,spikeAtWork(1,:),'b+',resampleBase,interpSpike,'r',resampleBase,interpSpike2,'k');
                %                 hold on
                %                 plot(1:21,spikeAtWork,'b+');
                %                 plot(1:21,spikeAtWork,'k--');
                %                 plot(2:20,centeredSpike,'r-');
                %                 hold off
                
            end % spikeTime
        end % el
    end % while ~isFinihed
    
end