function [projSpikes,eigenValues,eigenVectors,spikeTimes] = PCProj(dataPath, timeCommand, spikesTotal, covMatrix, averages, totSpikes)
    % Computes eigenvectors of electrode covariance matrices and computes spike projections
    
    
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
    
    nDims = projConfig.nDims;
    
    %% Creating data source
    dataSource = DataFileUpsampler(rawDataSource, meanTimeConstant, nLPoints, nRPoints);
    
    nElectrodes = dataSource.nElectrodes;
    disconnected = dataSource.disconnected;
    
    %% Setting up neighbor map
    % Subfunction encapsulates java use
    [adjacent,maxAdjacent] = catchAdjWJava( dataSource, projConfig.electrodeUsage);
    
    %% Data flow
    
    upSampRatio = dataSource.upSampleRatio;
    upSampStep = 1/upSampRatio;
    
    % interpolation bases
    resampleBase = 0:upSampStep:2;
    
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
            
            if dataSource.disconnected(el) || nSpikes == 0;
                continue
            end
            
            %% Check if buffer expansion is required
            if nSpikes > size(spikeBuffer,1);
                spikeBuffer = zeros(nSpikes,(nPoints-2) * maxAdjacent);
            end
            
            %% Store spike times
            spikeTimes{el}(currSpike(el):(currSpike(el) + nSpikes - 1)) = spikes{el};
            
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
                elAdj = adjacent{el}(elAdjIndex);
                spikeBuffer(1:nSpikes,((nPoints-2)*(elAdjIndex-1)+1):((nPoints-2)*elAdjIndex)) =...
                    reshape(dataSource.interpolant{elAdj}(interpPointsLin),s);
            end
            
            spikesTemp = spikeBuffer(1:nSpikes,1:(numel(adjacent{el})*(nPoints-2)));
            
            projSpikes{el}(currSpike(el):(currSpike(el) + nSpikes - 1),:) = ...
                bsxfun(@minus,spikesTemp,averages{el}) * eigenVectors{el};
            currSpike(el) = currSpike(el) + nSpikes;
            
            
            if false % Alignment debug plots
            %%
                clf
                for elAdjIndex = 1:numel(adjacent{el})
                    elAdj = adjacent{el}(elAdjIndex);
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
    end % while ~isFinihed
    
end