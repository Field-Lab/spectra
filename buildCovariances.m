function [covMatrix,averages,totSpikes] = buildCovariances(parameters, spikesTotal)
    % Build the covariance matrix for spikes around a given electrode
    % Input HashMap parameters should be the same given than for SpikeFindingM
    
    
    %% Imports
    import edu.ucsc.neurobiology.vision.electrodemap.*
    import edu.ucsc.neurobiology.vision.io.*
    import java.io.*
    
    %% Argument validation
    % Argument should be a java.util.HashMap<String,String> containing all relevant parameters for spike
    % finding
    validateattributes(parameters,{'java.util.HashMap'},{},'','parameters');
    p = parameters; % For concision
    
    
    %% Parsing and Storing Input HashMap
    % Note: Not so much a good input strategy
    % does not match so well CovarianceCalculator constructor
    
    rawDataSource = p.get('Raw_Data_Source'); % Actually at this point includes a command concatenated under the dataFileParser format: '.../data002(0-10)'
    sigmaPath = p.get('Sigma'); % .noise file
    outputPath = p.get('Analysis.Output_Path'); % Output path for the .spikes file
    
    meanTimeConstant = str2double(p.get('Mean Time Constant'));
    
    nLPoints = str2double(p.get('Analysis.Left Points'));
    nRPoints = str2double(p.get('Analysis.Right Points'));
    nPoints = nLPoints + nRPoints + 1;
    minError = str2double(p.get('Analysis.Minimization Error'));
    spikeUse = p.get('Analysis.Spike To Use');
    
    electrodeUsage = str2double(p.get('Analysis.Electrode Usage'));
    electrodeUsage = 1;
    
    %% Creating data source
    dataSource = DataFileUpsampler(rawDataSource, meanTimeConstant, nLPoints, nRPoints);
    
    
    %% Java electrodemap setup
    header = dataSource.rawDataFile.getHeader();
    packedArrayID = int32(header.getArrayID());
    
    electrodeMap = ElectrodeMapFactory.getElectrodeMap(packedArrayID);
    nElectrodes = electrodeMap.getNumberOfElectrodes();
    
    %% Setting up neighbor map
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
                
            end % spikeTime
        end % el
    end % while ~isFinished
    
    %% Normalize averages and covmatrix
    for el = 2:nElectrodes
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