function covMatrix = buildCovariances(parameters, spikeFileName)
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
    
    bufferLengthInSamples = samplingRate;
    
    while ~isFinished % stopSample should be the first sample not loaded
        bufferStart = max(startSample, lastSampleLoaded + 1 - nLPoints); % Buffer beginning (inclusive)
        bufferEnd = min(lastSampleLoaded + bufferLengthInSamples + nRPoints + 1, stopSample); % Buffer end (exclusive)
        if bufferEnd == stopSample
            isFinished = true;
        end
        lastSampleLoaded = bufferEnd - nRPoints - 1;
        
        rawData = rawDataFile.getData(bufferStart, bufferEnd - bufferStart)';
        
        spikes = spikeFile.getSpikesTimesUntil(bufferStart + nLPoints, bufferEnd - nRPoints);
        for el = 2:nElectrodes % Could parallel here, but actually slower due to the IO cost of sending to each worker.
            % The better part would be to change the while loop to a smarter parfor loop.
            for spikeTime = spikes{el}'
                spikeAtWork = ...
                    rawData(adjacent{el}+1,(spikeTime-nLPoints-bufferStart+1):(spikeTime+nRPoints-bufferStart+1));
%                     plot(spikeAtWork');
%                     pause(0.01);
%                     input(['el: ',num2str(el),' ; spikeTime: ',num2str(spikeTime)]);
            end
        end
    end
    
    
end