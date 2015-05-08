function [covMatrix,averages,totSpikes] = buildCovariances(parameters, spikeFileName)
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
    
    %% Creating data source
    % We are not listening to any sampleInputStream
    % --> Setting up a data source straight from file
    parser = DataFileStringParser(rawDataSource);
    datasets = parser.getDatasets();
    rawDataFile = RawDataFile(File(char(datasets(1))));
    startTimes = parser.getStartTimes();
    stopTimes = parser.getStopTimes();
    
    header = rawDataFile.getHeader();
    samplingRate = header.getSamplingFrequency();
    
    startSample = startTimes(1) * samplingRate;
    stopSample = stopTimes(1) * samplingRate;
    
    totalSamples = stopSample - startSample;
    
    alpha = 1 / (meanTimeConstant * samplingRate);
    
    %% Creating spike source
    spikeFile = SpikeFile(spikeFileName);
    
    %% Java electrodemap setup
    packedArrayID = int32(header.getArrayID());
    binPackedArrayID = dec2bin(packedArrayID,32);
    arrayID = bin2dec(binPackedArrayID(17:32));
    arrayPart = bin2dec(binPackedArrayID(9:16));
    arrayNParts = bin2dec(binPackedArrayID(1:8));
    if arrayPart == 0
        arrayPart = 1;
        arrayNParts = 1;
    end
    
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
    
    %% Data strategy:
    % Load buffers of 1 sec for spikes
    % Load buffers of 1 sec + edge points for matrix
    
    lastSampleLoaded = startSample-1;
    isFinished = false;
    
    bufferLengthInSamples = 2500; %samplingRate;
    validateattributes(bufferLengthInSamples,{'numeric'},{'scalar','integer','>',0},'','bufferLengthInSamples');
    
    % Initialize filter state for data filtering - 1st order highpass IIR
    filterState = zeros(1,513);
    bFilter = (1-alpha)*[1,-1];
    aFilter = [1,alpha-1];
    
    % Assigning a UniformSpline java object -- java hybrid debug purposes
    % splineObject = edu.ucsc.neurobiology.vision.math.UniformSpline(nPoints);
    
    % interpolation bases
    interpBase = 1:nPoints;
    resampleBase = nLPoints:0.01:(nLPoints+2);
    
    % Initializing covariance matrices
    covMatrix = cell(nElectrodes,1);
    averages = cell(nElectrodes,1);
    totSpikes = zeros(nElectrodes,1);
    
    for i = 2:nElectrodes
        covMatrix{i} = zeros((nPoints-2) * numel(adjacent{i}));
        averages{i} = zeros(1,(nPoints - 2) * numel(adjacent{i}));
    end
    
    %%%
    stopSample = 20000;
    upSampRatio = 100;
    upSampStep = 1/upSampRatio;
    %%%
    while ~isFinished % stopSample should be the first sample not loaded
        %% Load samples
        bufferStart = max(startSample, lastSampleLoaded + 1 - nLPoints) % Buffer beginning (inclusive)
        bufferEnd = min(lastSampleLoaded + bufferLengthInSamples + nRPoints + 1, stopSample) % Buffer end (exclusive)
        if bufferEnd == stopSample
            isFinished = true;
        end
        lastSampleLoaded = bufferEnd - nRPoints - 1;
        
        % Filter
        [rawData, filterState] = filter(bFilter, aFilter, ...
            single(rawDataFile.getData(bufferStart, bufferEnd - bufferStart)'),...
            filterState, 2);
        
        %% per electrode interp model. Too many losses in ppval otherwise.
        upSampPoints = 1:upSampStep:size(rawData,2);
        upSampData = [];
        upSampData(nElectrodes,numel(upSampPoints)) = single(0);
        
        for el = nElectrodes:-1:1
            interpModel{el} = griddedInterpolant(1:size(rawData,2),rawData(el,:),'spline');
            upSampData(el,:) = interpModel{el}(upSampPoints);
        end
        
        %% Load Spikes
        spikes = spikeFile.getSpikesTimesUntil(bufferStart + nLPoints, bufferEnd - nRPoints);
        % If no spikes at all are loaded, skip iteration
        % Required as by Matlab cast spikes is empty 513x0 and not a cell array in that case
        if size(spikes,2) == 0
            continue
        end
        
        %% Process by electrodes
        % Could parallel here, but actually slower due to the IO cost of sending to each worker.
        % The better part would be to change the while loop to a smarter parfor loop.
        for el = 2:nElectrodes
            %% Process each spike
            for spikeTime = spikes{el}'
                % Load master spike
                interpSpike = upSampData(el,round(upSampRatio*(resampleBase + double(spikeTime) - bufferStart - nLPoints - 1))+1);
                
                % Find minimum and compute associated resample points
                offset = (find(interpSpike == min(interpSpike),1)-1)/100;
                interpPoints = (1:(nPoints-2)) + offset;
                
                % Load realigned spikes, master + neighbors
                centeredSpike = upSampData(adjacent{el}+1,round(upSampRatio*(interpPoints + double(spikeTime) - bufferStart - nLPoints - 1))+1)';
                
                % Update info
                totSpikes(el) = totSpikes(el) + 1;
                averages{el} = averages{el} + centeredSpike(:)';
                covMatrix{el} = covMatrix{el} + centeredSpike(:) * centeredSpike(:)';
                
%                 plot(interpBase,spikeAtWork(1,:),'b+',resampleBase,interpSpike,'r');
%                 hold on
%                 plot(1:21,spikeAtWork,'b+');
%                 plot(1:21,spikeAtWork,'k--');
%                 plot(2:20,centeredSpike,'r-');
%                 hold off
            end
        end
    end
    
    %% Normalize averages and covmatrix
    for el = 2:nElectrodes
        if totSpikes(el) >= 2
            covMatrix{el} = (covMatrix{el} - averages{el}' * averages{el} / totSpikes(el))/(totSpikes(el)-1);
            %             covMatrix{el} = (covMatrix{el} / totSpikes(el))/(totSpikes(el)-1);
            averages{el} = averages{el}/totSpikes(el);
        end
    end
    
end