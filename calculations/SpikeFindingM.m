function [spikes,ttlTimes,nSamples] = SpikeFindingM(dataPath, saveFolder, timeCommand, sigmaPath )
    %SPIKEFINDINGM Finds spikes over a dataset, using a SpikeFinder instance
    %   for thresholdind and time-checking operations
    %
    % Inputs:
    %   dataPath: path to the raw data files
    %   saveFolder: path to the analysis folder
    %   timeCommand: if any, time tag of the data subset to analyze
    %   sigmaPath: path to the data noise floor values file, in ascii format
    %
    % Outputs:
    %   spikes: nSpikes x 3 array containing all spikes occurences on the dataset
    %       1st col (sorted) spike times - 2nd col electrodes - 3rd col max amplitudes
    %   ttlTimes: nTtl x 1 array containing all ttlTimes for the dataset
    %
    % Author -- Vincent Deo -- Stanford University -- August 27, 2015
    
    %% Setup
    % Argument validation
    if ~(exist(dataPath,'file') == 2 || exist(dataPath,'file') == 7)
        throw(MException('','SpikeFinding: data folder|file does not exist'));
    end
    if ~(exist(saveFolder,'file') == 7)
        throw(MException('','SpikeFinding: output folder does not exist'));
    end
    if ~(exist(sigmaPath,'file') == 2)
        throw(MException('','SpikeFinding: noise file does not exist'));
    end
    
    % Load spike finding configuration
    config = mVisionConfig();
    spikeConfig = config.getSpikeConfig();
    
    % Parameters and datasource setup
    rawDataSource = [dataPath,timeCommand];
    
    spikeThreshold = spikeConfig.spikeThreshold;
    ttlThreshold = spikeConfig.ttlThreshold;
    meanTimeConstant = spikeConfig.meanTimeConstant;
    
    dataSource = DataFileUpsampler(rawDataSource, meanTimeConstant, 0, 0);
    nElectrodes = dataSource.nElectrodes;
   
    % Load noise values
    sigma = load(sigmaPath,'-ascii');
    validateattributes(sigma,{'numeric'},{'ncols',1,'nrows',nElectrodes},'','loaded sigma file');

    sigma(1) = spikeConfig.defaultTtlSigma; % TTL electrode
    sigma = sigma * spikeThreshold; % Apply noise floor multiplier
     
    % Instantiate Spiker Finder
    spikeFinderM = SpikeFinderM(sigma, ttlThreshold, dataSource);
    
    %% Spike Finding
    
    % Initialization
    [s,f] = dataSource.loadNextBuffer(dataSource.samplingRate,false); % Get a read of 1 second
    if spikeConfig.debug
        disp(['Buffer ',num2str(s),' - ',num2str(f)]);
    end
            
    dataSource.filterState = -spikeFinderM.initialize()'; % Get filter initializing values - av. of first second
    
    if size(dataSource.rawData,2) ~= dataSource.samplingRate
        throw(MException('',...
            'SpikeFindingM:Total number of samples is insufficient to initialize SpikeFinder'));
    end
    
    firstIter = true;
    
    spikes = zeros(0,3);
    
    % Main Buffer loop
    while ~dataSource.isFinished % stopSample should be the first sample not loaded
        if ~firstIter % Load a new buffer - filtering included
            [s,f] = dataSource.loadNextBuffer();
            if spikeConfig.debug
                disp(['Buffer ',num2str(s),' - ',num2str(f)]);
            end
        else % Forcefiltering of unfiltered already loaded initialization buffer
            firstIter = false;
            dataSource.forceFilter(true);
        end
        
        % Work the spike finder
        spikes = [spikes;spikeFinderM.processBuffer()];
    end % buffer loop
    
    % Sort output by time, isolate TTLs
    spikes = sortrows(spikes,1);
    ttlTimes = spikes(spikes(:,2) == 1,1);
    nSamples = dataSource.stopSample - dataSource.startSample;
end