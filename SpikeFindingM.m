function [spikes,ttlTimes] = SpikeFindingM(dataPath, saveFolder, timeCommand, sigmaPath )
    %SPIKEFINDINGM Matlab implementation of the Spike Finding algorithm
    %   Takes on after noise finding has been done (or not)
    % Computes the spikesproperties over the electrode array and stores
    % them in a .spikes file
    % In a format matching the current behavior of vision.
    
    import edu.ucsc.neurobiology.vision.electrodemap.*
    import edu.ucsc.neurobiology.vision.io.*
    import java.io.*
    
    %% Argument validation
    if ~(exist(dataPath,'file') == 2 || exist(dataPath,'file') == 7)
        throw(MException('','SpikeFinding: data folder|file does not exist'));
    end
    if ~(exist(saveFolder,'file') == 7)
        throw(MException('','SpikeFinding: output folder does not exist'));
    end
    if ~(exist(sigmaPath,'file') == 2)
        throw(MException('','SpikeFinding: noise file does not exist'));
    end
    
    %% Loading spike finding configuration
    config = mVisionConfig();
    spikeConfig = config.getSpikeConfig();
    
    %% Parsing and Storing input HashMap
    rawDataSource = [dataPath,timeCommand];
    
    spikeThreshold = spikeConfig.spikeThreshold;
    ttlThreshold = spikeConfig.ttlThreshold;
    meanTimeConstant = spikeConfig.meanTimeConstant;
    
    %% Creating data source
    dataSource = DataFileUpsampler(rawDataSource, meanTimeConstant, 0, 0);
    
    %% Java electrodemap setup
%     header = dataSource.rawDataFile.getHeader();
%     packedArrayID = int32(header.getArrayID());
    
%     electrodeMap = ElectrodeMapFactory.getElectrodeMap(packedArrayID);
    nElectrodes = dataSource.nElectrodes;
   
    %% Get sigmas
    sigma = load(sigmaPath,'-ascii');
    validateattributes(sigma,{'numeric'},{'ncols',1,'nrows',nElectrodes},'','loaded sigma file');

    sigma(1) = spikeConfig.defaultTtlSigma; % TTL electrode
    sigma = sigma * spikeThreshold;
     
    %% Create the Spiker Finder and heirs
    spikeFinderM = SpikeFinderM(sigma, ttlThreshold, meanTimeConstant, dataSource);
    
    %% Spike Finding
    
    % Initialization
    [s,f] = dataSource.loadNextBuffer(dataSource.samplingRate,false); % Get a read of 1 second
    if spikeConfig.debug
        disp(['Buffer ',num2str(s),' - ',num2str(f)]);
    end
            
    dataSource.filterState = -spikeFinderM.initialize()';
    
    if size(dataSource.rawData,2) ~= dataSource.samplingRate
        throw(MException('',...
            'SpikeFindingM:Total number of samples is insufficient to initialize SpikeFinder'));
    end
    
    firstIter = true;
    
    spikes = zeros(0,3);
    
    while ~dataSource.isFinished % stopSample should be the first sample not loaded
        if ~firstIter
            [s,f] = dataSource.loadNextBuffer();
            if spikeConfig.debug
                disp(['Buffer ',num2str(s),' - ',num2str(f)]);
            end
        else
            firstIter = false;
            dataSource.forceFilter(true);
        end
        
        spikes = [spikes;spikeFinderM.processBuffer()];
    end
    
    spikes = sortrows(spikes,1);
    ttlTimes = spikes(spikes(:,2) == 1,1);
    
end