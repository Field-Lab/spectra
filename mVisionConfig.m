classdef mVisionConfig
    %MVISIONCONFIG Configuration class for mvision scripts
    %   Contains parameters to run ANF pipeline
    %   Without using a java/xml vision config object/file
    %
    % Author -- Vincent Deo -- Stanford University -- August 7, 2015
    
    properties (SetAccess = immutable, GetAccess = private)
        % General
        debug = false;
        
        % Parallel (multi-file at parallelCaller level, clustering loop, ...)
        nWorkers                = 32 % 01/04/16: Poorly managed, pool times out before reaching clustering
        
        %% Data Source management
        bufferMaxSize           = 16384 % samples
        upSampleRatio           = 16    % samples
        sampleRate              = 20000 % Hz
        
        %% Raw Data Noise Evaluation
        noiseTime               = 5     % seconds
        noiseTimeToSkip         = 5     % seconds
        
        %% Spike Finding properties
        defaultTtlSigma         = 100   % Amplitude
        spikeThreshold          = 4     % Amplitude factor
        ttlThreshold            = 1000  % Amplitude factor
        meanTimeConstant        = 0.01  % Seconds
        
        minSpikeSeparation      = 5.0   % samples
        maxSpikeWidth           = 50.0  % samples
        
        %% Covariance Calculation properties
        nLPoints                = 5     % samples
        nRPoints                = 15    % samples
        % 0 - 1 electrode || 1 - 7 electrodes || 2 - 19 electrodes
        electrodeUsage          = 1
        
        whitening               = true  % whiten relative to noise floor on electrodes
        expectedNoiseEvents     = 10000;
        minNoiseSpacing         = 20;
        
        covSpikeBufferSize      = 100;
        
        %% Projections properties
        projNDimensions         = 5
        
        %% Clustering management properties
        maxClusteringTries      = 10    %
        minSpikesForClustering  = 100   %
        
        %% Spectral clustering properties
        spikesForLaplacian      = 2500  %
        spikesForKMeans         = 20000 %
        nthNeighbor             = 5     %
        maxEigenVectAndClusters = 15    % Maximum dimension of the eigenvector subspace in which we cluster - Actual dimension is numClusters
        eigsMaxIter             = 300   %
        eigsTol                 = eps   %
        qualityTol              = 0.015 %
        kMeansReplicas          = 5     %
        kMeansMaxIterations     = 200   %
        reprojBufferSize        = 1000  %
        maxOutlierFraction      = 0.01  %
        
        %% Neuron cleaning properties
        minSpikes               = 100   % spikes
        maxContamination        = 0.1   %
        
        EILeftPoints            = 20    % samples
        EIRightpoints           = 60    % samples
        EISpikesToAverage       = 3000  % spikes
        EInThreads              = 4     % threads
        
        EIMergeThresholdWithin  = 0.1   % Normalized metric
        EIMergeThresholdAcross  = 0.1   % Normalized metric
        
        EIGlobalMinWindow       =[18 44]% Sample window in which perform global min EI realignment
        
    end % properties
    
    methods
        
        function dataConfig = getDataConfig(obj)
            dataConfig.bufferMaxSize = obj.bufferMaxSize;
            dataConfig.upSampleRatio = obj.upSampleRatio;
            dataConfig.sampleRate = obj.sampleRate;
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
            
            covConfig.whitening = obj.whitening;
            covConfig.noiseEvents = obj.expectedNoiseEvents;
            covConfig.noiseSpacing = obj.minNoiseSpacing;
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
            clustConfig.minSpikes = obj.minSpikesForClustering;
            clustConfig.maxTries = obj.maxClusteringTries;
        end % getClustConfig
        
        function specConfig = getSpectralConfig(obj)
            specConfig.debug = obj.debug;
            
            specConfig.nthNeighbor = obj.nthNeighbor;
            specConfig.spikesForLaplacian = obj.spikesForLaplacian;
            specConfig.spikesForKMeans = obj.spikesForKMeans; 
            specConfig.maxEV = obj.maxEigenVectAndClusters;
            specConfig.qualityTol = obj.qualityTol;
            specConfig.maxOutlierFraction = obj.maxOutlierFraction;
            
            specConfig.eigsMaxIter = obj.eigsMaxIter;
            specConfig.eigsTol = obj.eigsTol;
            
            specConfig.reprojBufferSize = obj.reprojBufferSize;
            specConfig.kMeansReplicas = obj.kMeansReplicas;
            specConfig.kMeansIterations = obj.kMeansMaxIterations;
        end % getSpectralConfig
        
        function parConfig = getParConfig(obj)
            parConfig.nWorkers= obj.nWorkers;
        end % getParconfig
        
        function cleanConfig = getCleanConfig(obj)
            cleanConfig.minSpikes = obj.minSpikes;
            cleanConfig.maxCont = obj.maxContamination;
            
            cleanConfig.EITC = obj.meanTimeConstant;
            cleanConfig.EILP = obj.EILeftPoints;
            cleanConfig.EIRP = obj.EIRightpoints;
            cleanConfig.EISp = obj.EISpikesToAverage;
            cleanConfig.EInThreads = obj.EInThreads;
            
            cleanConfig.eiThrW = obj.EIMergeThresholdWithin;
            cleanConfig.eiThrA = obj.EIMergeThresholdAcross;
            
            cleanConfig.globMinWin = obj.EIGlobalMinWindow;
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

