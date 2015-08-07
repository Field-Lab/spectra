function [projSpikes,eigenValues,eigenVectors,spikeTimes] = PCProj(dataPath, timeCommand, spikesTotal, covMatrix, averages, totSpikes)
    % Computes eigenvectors of electrode covariance matrices and computes spike projections
    
    
    %% Argument validation
    if ~(exist(dataPath,'file') == 2 || exist(dataPath,'file') == 7)
        throw(MException('','ProjCalculations: data folder|file does not exist'));
    end
    
    %% Load projections configuration
    config = mVisionConfig();
    projConfig = config.getProjConfig();
    
    rawDataSource = [dataPath,timeCommand];
    
    meanTimeConstant = projConfig.meanTimeConstant;
    
    nLPoints = projConfig.nLPoints;
    nRPoints = projConfig.nRPoints;
    
    nDims = projConfig.nDims;
    
    %% Creating data source
    dataSource = DataFileUpsampler(rawDataSource, meanTimeConstant, nLPoints, nRPoints);
    
    nElectrodes = dataSource.nElectrodes;
    disconnected = dataSource.disconnected;
    
    aligner = spikeAligner(dataSource);
    
    %% Projections storage and init
    eigenValues = cell(nElectrodes,1);
    eigenVectors = cell(nElectrodes,1);
    projSpikes = cell(nElectrodes,1);
    spikeTimes = cell(nElectrodes,1);
    currSpike = ones(nElectrodes,1);
    
    whitener = cell(nElectrodes,1);
    
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
        
        whitener{el} = fliplr(v) * diag(e.^-0.5) * fliplr(v)';
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
            
            if dataSource.disconnected(el) || nSpikes == 0;
                continue
            end

            %% Store spike times
            spikeTimes{el}(currSpike(el):(currSpike(el) + nSpikes - 1)) = spikes{el};
            
            %% Align all spikes with the aligner            
            spikesTemp = aligner.alignSpikes(el,spikes{el});
            
            % Standard
            projSpikes{el}(currSpike(el):(currSpike(el) + nSpikes - 1),:) = ...
                bsxfun(@minus,spikesTemp,averages{el}) * eigenVectors{el};
            
%             % Whitened
%             spikesTemp2 = sqrt(nSpikes) * bsxfun(@minus,spikesTemp,averages{el}) * whitener{el};
%             projSpikes{el}(currSpike(el):(currSpike(el) + nSpikes - 1),:) = spikesTemp2(:,1:nDims);
                
            currSpike(el) = currSpike(el) + nSpikes;
            
            
        end % el
    end % while ~isFinihed
    
end