classdef mVisionConfig
    %MVISIONCONFIG Configuration class for mvision scripts
    %   Contains parameters to run ANF pipeline
    %   Without using a java/xml vision config object/file
    
    
    properties (SetAccess = immutable, GetAccess = public)
        % General
        debug = false;
        
        %% Data Source management
        bufferMaxSize           = 4096  % samples
        upSampleRatio           = 16    % samples
        
        %% Raw Data Noise Evaluation
        noiseTime               = 5     % seconds
        noiseTimeToSkip         = 5     % seconds
        
        %% Spike Finding properties
        defaultTtlSigma         = 100   % Amplitude
        spikeThreshold          = 3     % Amplitude factor
        ttlThreshold            = 1000  % Amplitude factor
        meanTimeConstant        = 0.01  % Seconds
        
        %% Covariance Calculation properties
        nLPoints                = 5     % samples
        nRPoints                = 15    % samples
        % 0 - 1 electrode || 1 - 7 electrodes || 2 - 19 electrodes
        electrodeUsage          = 1
        
        covSpikeBufferSize      = 100;
        
        %% Projections properties
        projNDimensions         = 3
        
        %% Clustering properties
        opticsSubsetMaxSize     = 5000  % PC space points
        binsPerDimension        = 30    % N-dimensional bins
        opticsDensityFactor     = 5     %
        overlapFactorForDiscard = 2     % Distance between means in sigmas
        maxGaussians            = 8     %
        maxEMIter               = 300   %
        regularizationValue     = 0.01  %
        
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
            clustConfig.opticsSubsetMaxSize = obj.opticsSubsetMaxSize;
            clustConfig.binsPerDimension = obj.binsPerDimension;
            clustConfig.opticsDensityFactor = obj.opticsDensityFactor;
            clustConfig.overlapFactorForDiscard = obj.overlapFactorForDiscard;
            clustConfig.maxGaussians = obj.maxGaussians;
            clustConfig.maxEMIter = obj.maxEMIter;
            clustConfig.regVal = obj.regularizationValue;
        end % getClustConfig
        
    end % methods
    
end

