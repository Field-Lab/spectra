classdef mVisionConfig
    %MVISIONCONFIG Configuration class for mvision scripts
    %   Contains parameters to run ANF pipeline
    %   Without using a java/xml vision config object/file
    
    
    properties (SetAccess = immutable, GetAccess = private)
        % General
        debug = true;
        
        % Parallel (multi-file at parallelCaller level)
        nWorkers                = 4
        
        %% Data Source management
        bufferMaxSize           = 65536 % samples
        upSampleRatio           = 16    % samples
        
        %% Raw Data Noise Evaluation
        noiseTime               = 5     % seconds
        noiseTimeToSkip         = 5     % seconds
        
        %% Spike Finding properties
        defaultTtlSigma         = 100   % Amplitude
        spikeThreshold          = 3     % Amplitude factor
        ttlThreshold            = 1000  % Amplitude factor
        meanTimeConstant        = 0.01  % Seconds
        
        minSpikeSeparation      = 5.0   % samples
        maxSpikeWidth           = 50.0  % samples
        
        %% Covariance Calculation properties
        nLPoints                = 5     % samples
        nRPoints                = 15    % samples
        % 0 - 1 electrode || 1 - 7 electrodes || 2 - 19 electrodes
        electrodeUsage          = 1
        
        covSpikeBufferSize      = 100;
        
        %% Projections properties
        projNDimensions         = 5
        
        %% Clustering properties
        maxGaussians            = 8     %
        maxEMIter               = 300   %
        regularizationValue     = 0.001 %
        belongProbability       = 0.6   %
        
        % Neuron cleaning properties
        minSpikes               = 100   % spikes
        maxContamination        = 0.1   %
        coincidenceTime         = 10    % samples
        maxCorrelation          = 0.25  %
        
        
    end % properties
    
    methods
        
        function dataConfig = getDataConfig(obj)
            dataConfig.bufferMaxSize = obj.bufferMaxSize;
            dataConfig.upSampleRatio = obj.upSampleRatio;
        end % dataConfig
        
        function noiseConfig = getNoiseConfig(obj)
            noiseConfig.debug = obj.debug;
            
            noiseConfig.time = obj.noiseTime;
            noiseConfig.timeToSkip = obj.noiseTimeToSkip;
        end % getNoiseConfig
        
        function spikeConfig = getSpikeConfig(obj)
            spikeConfig.debug = obj.debug;
            
            spikeConfig.spikeThreshold = obj.spikeThreshold;
            spikeConfig.ttlThreshold = obj.ttlThreshold;
            spikeConfig.meanTimeConstant = obj.meanTimeConstant;
            spikeConfig.defaultTtlSigma = obj.defaultTtlSigma;
            
            spikeConfig.minSpikeSeparation = obj.minSpikeSeparation;
            spikeConfig.maxSpikeWidth = obj.maxSpikeWidth;
        end % getSpikeConfig
    
        function covConfig = getCovConfig(obj)
            covConfig.debug = obj.debug;
            
            covConfig.meanTimeConstant = obj.meanTimeConstant;
            covConfig.nLPoints = obj.nLPoints;
            covConfig.nRPoints = obj.nRPoints;
            covConfig.electrodeUsage = obj.electrodeUsage;
            covConfig.spikeBufferSize = obj.covSpikeBufferSize;
        end % getCovConfig
        
        function projConfig = getProjConfig(obj)
            projConfig.debug = obj.debug;
            
            projConfig.meanTimeConstant = obj.meanTimeConstant;
            projConfig.nLPoints = obj.nLPoints;
            projConfig.nRPoints = obj.nRPoints;
            projConfig.electrodeUsage = obj.electrodeUsage;
            projConfig.nDims = obj.projNDimensions;
        end % getProjConfig
        
        function clustConfig = getClustConfig(obj)
            clustConfig.debug = obj.debug;
            
            clustConfig.nDims = obj.projNDimensions;
            clustConfig.maxGaussians = obj.maxGaussians;
            clustConfig.maxEMIter = obj.maxEMIter;
            clustConfig.regVal = obj.regularizationValue;
            clustConfig.clusterProb = obj.belongProbability;
        end % getClustConfig
        
        function parConfig = getParConfig(obj)
            parConfig.nWorkers= obj.nWorkers;
        end % getParconfig
        
        function cleanConfig = getCleanConfig(obj)
            cleanConfig.minSpikes = obj.minSpikes;
            cleanConfig.maxCont = obj.maxContamination;
            cleanConfig.coincTime = obj.coincidenceTime;
            cleanConfig.maxCorr = obj.maxCorrelation;
        end
        
        function pools = buildElPools(obj,nElectrodes)
            cuts = round((0:obj.nWorkers) .* nElectrodes ./ obj.nWorkers) + 1;
            pools = cell(obj.nWorkers,1);
            for p = 1:obj.nWorkers
               pools{p} = cuts(p):(cuts(p+1)-1); 
            end
        end
        
    end % methods
    
end

