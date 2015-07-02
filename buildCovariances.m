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
    resampleBase = 0:upSampStep:2;
    
    % Initializing covariance matrices
    covMatrix = cell(nElectrodes,1);
    averages = cell(nElectrodes,1);
    totSpikes = zeros(nElectrodes,1);
        
    for i = 2:nElectrodes
        covMatrix{i} = zeros((nPoints-2) * numel(adjacent{i}));
        averages{i} = zeros(1,(nPoints - 2) * numel(adjacent{i}));
    end
    
    spikeBuffer = zeros(1,(nPoints-2) * maxAdjacent);
    
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
            
            if dataSource.disconnected(el) || nSpikes == 0
                continue
            end            
            
            %% Check if buffer expansion is required
            if nSpikes > size(spikeBuffer,1);
                spikeBuffer = zeros(nSpikes,(nPoints-2) * maxAdjacent);
            end
            
            
            %% Process all spikes
            interpIndex = bsxfun(@plus,resampleBase,double(spikes{el}) -  bufferStart);
            interpSpikes = dataSource.interpolant{el}(interpIndex(:));
            interpSpikes = reshape(interpSpikes,size(interpIndex));
            
            [~,offset] = min(interpSpikes,[],2);
            interpPoints = bsxfun(@plus,...
                (offset-1)/upSampRatio + double(spikes{el}) - bufferStart - nLPoints,...
                1:(nPoints-2));
            interpPointsLin = interpPoints(:);
            
            s = size(interpPoints);
            for elAdjIndex = 1:numel(adjacent{el})
                elAdj = adjacent{el}(elAdjIndex) + 1;
                spikeBuffer(1:nSpikes,((nPoints-2)*(elAdjIndex-1)+1):((nPoints-2)*elAdjIndex)) =...
                    reshape(dataSource.interpolant{elAdj}(interpPointsLin),s);
            end
            
            totSpikes(el) = totSpikes(el) + nSpikes;
            
            spikesTemp = spikeBuffer(1:nSpikes,1:(numel(adjacent{el})*(nPoints-2)));
            
            averages{el} = averages{el} + sum(spikesTemp,1);
            covMatrix{el} = covMatrix{el} + spikesTemp' * spikesTemp;
            
            if false % Alignment debug plots
            %%
                clf
                for elAdjIndex = 1:numel(adjacent{el})
                    elAdj = adjacent{el}(elAdjIndex) + 1;
                    hold on
                    plot(1:size(dataSource.rawData,2),dataSource.rawData(elAdj,:)+(elAdjIndex-1)*150,'k+');
                    plot(1:0.05:size(dataSource.rawData,2),...
                        dataSource.interpolant{elAdj}(1:0.05:size(dataSource.rawData,2))+(elAdjIndex-1)*150,'b--')
                    for sp = 1:nSpikes
                        plot(interpPoints(sp,:),...
                        spikesTemp(sp,((elAdjIndex-1)*(nPoints-2) + 1):(elAdjIndex*(nPoints-2)))+(elAdjIndex-1)*150,'r+-');
                    end
                    offset
                    hold off
                end
            end
            
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