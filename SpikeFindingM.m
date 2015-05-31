function SpikeFindingM( parameters )
    %SPIKEFINDINGM Matlab implementation of the Spike Finding algorithm
    %   Takes on after noise finding has been done (or not)
    % Computes the spikesproperties over the electrode array and stores
    % them in a .spike file
    % In a format matching the current behavior of vision.
    
    import edu.ucsc.neurobiology.vision.electrodemap.*
    import edu.ucsc.neurobiology.vision.io.*
    import java.io.*
    
    %% Argument validation
    % Argument should be a java.util.HashMap<String,String> containing all relevant parameters for spike
    % finding
    validateattributes(parameters,{'java.util.HashMap'},{},'','parameters');
    p = parameters; % For concision
    
    %% Parsing and Storing input HashMap
    rawDataSource = p.get('Raw_Data_Source'); % Actually at this point includes a command concatenated under the dataFileParser format: '.../data002(0-10)'
    sigmaPath = p.get('Sigma'); % .noise file
    outputPath = p.get('Analysis.Output_Path'); % Output path for the .spikes file
    
    spikeThreshold = str2double(p.get('Spike Threshold'));
    ttlThreshold = str2double(p.get('TTL Threshold'));
    meanTimeConstant = str2double(p.get('Mean Time Constant'));
    
    %% Creating data source
    dataSource = DataFileUpsampler(rawDataSource, meanTimeConstant, 0, 0);
    
    %% Java electrodemap setup
    header = dataSource.rawDataFile.getHeader();
    packedArrayID = int32(header.getArrayID());
    
    electrodeMap = ElectrodeMapFactory.getElectrodeMap(packedArrayID);
    nElectrodes = electrodeMap.getNumberOfElectrodes();
    
    disconnected = electrodeMap.getDisconnectedElectrodesList();
    
    %% Get sigmas
    % sigma = getSigmas(sigmaPath, nElectrodes)
    % Implemeting functionality of SpikeFinding.getSigmas(String fileNameOrValue, nElectrodes);
    % TODO Check better and throw an error
    if isempty(str2num(sigmaPath)) % sigmaPath is a path string
        sigma = load(sigmaPath,'-ascii');
        validateattributes(sigma,{'numeric'},{'ncols',1,'nrows',nElectrodes},'','loaded sigma file');
        % File path error
    else % sigmaPath is a value
        sigma = sigmaPath * ones(nElectrodes,1);
    end
    sigma(1) = 100; % TTL electrode
    sigma = sigma * spikeThreshold;
    
    
    %% Raw Data Saving
    % Skipped - this version meant to work on already recorded data
    
    %% Create the Spiker Finder and heirs
    spikeFinderM = SpikeFinderM(electrodeMap, sigma, ttlThreshold, meanTimeConstant);
    spikeBufferM = SpikeBufferM();
    spikeSaverM  = SpikeSaverM(header, outputPath, meanTimeConstant, spikeThreshold);
    
    %% Spike Finding
    
    % Initialization
    [s,f] = dataSource.loadNextBuffer(dataSource.samplingRate); % Get a read of 1 second
    % lastSampleLoaded = samplingRate; % Do not update - after initializing resend all
    % initialization samples, with updated means.
    spikeFinderM.initialize(dataSource.rawData);
    
    if size(dataSource.rawData,2) ~= dataSource.samplingRate
        throw(MException('SpikeFindingM',...
            'Total number of samples is insufficient to initialize SpikeFinder'));
    end
    
    firstIter = true;
    currentSample = 0;
    
    while ~dataSource.isFinished % stopSample should be the first sample not loaded
        if ~firstIter
            [s,f] = dataSource.loadNextBuffer();
        else
            firstIter = false;
        end
        
        for i = 1:size(dataSource.rawData,2)
            currentSample = currentSample + 1;
            spikeBufferM.addSpikes(spikeFinderM.processSample(double(dataSource.rawData(:,i))));
            s = spikeBufferM.getSpikes(currentSample);
            spikeSaverM.processSpikes(s);
        end
    end
    
    spikeSaverM.processSpikes(spikeBufferM.getAllSpikes());
    spikeSaverM.finishSpikeProcessing();
    
end