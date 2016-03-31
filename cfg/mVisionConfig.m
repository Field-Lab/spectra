classdef mVisionConfig < handle
    %MVISIONCONFIG Configuration class for mvision scripts
    %   Contains parameters to run ANF pipeline
    %   Without using a java/xml vision config object/file
    %
    % Author -- Vincent Deo -- Stanford University -- August 7, 2015
    
    properties (SetAccess = immutable, GetAccess = public)
        %% General
        debug = false;
        
        %% Parallel (multi-file at parallelCaller level, clustering loop, ...)
        nWorkers                = 12
        
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
        expectedNoiseEvents     = 10000
        minNoiseSpacing         = 20
        
        covSpikeBufferSize      = 100
        
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
        
        %% Raw computations properties
        computeRawEIs           = true
        computeRawSTAs          = true
        keepNeuronsRaw          = false
        
        computeFinalEIs         = true
        computeFinalSTAs        = true
        computeFinalParams      = true
        
        %% Neuron cleaning properties
        minSpikes               = 100   % spikes
        maxContamination        = 0.1   %
        
        EILeftPoints            = 20    % samples
        EIRightPoints           = 60    % samples
        EISpikesToAverage       = 3000  % spikes
        EInThreads              = 4     % threads
        
        EIResamplePitch         = 0.1  % samples
        EIMergeThresholdWithin  = 0.1  % Normalized metric
        EIMergeThresholdGlobal  = 0.1  % Normalized metric
        
        EIGlobalMinWindow       = [18 44]% Sample window in which perform global min EI realignment
        
        STADepth                = 30   % Frames
        STAOffset               = -1   % Frames
        STAnSpikes              = 1000000000 % Spikes
        STACalcAtOnce           = 1000 % Spikes
        STAnThreads             = 10   % Threads
        
        paramsnThreads          = 10   % Threads
        
    end % properties
    
    methods
        
        % Constructor
        % Takes in a string that is the config path tag in STATIC_CONFIG_LIST
        %
        % The constructor parses the config list file
        % Then calls the destination file as a m source file
        % and will add 'obj.' at the beginning of each line
        %
        % A custom configuration file is a list of instructions (and comments)
        % of the form _parameter_ =  value;
        % Where parameter is a valid property name from the list above
        function obj = mVisionConfig(varargin)
            narginchk(0,1);
            if nargin == 0 || ( ischar(varargin{1}) && numel(varargin{1}) == 0) % case ''
                return;
            end
            
            cfgListPath = './cfg/STATIC_CONFIG_LIST'; % Path to STATIC_CONFIG_LIST from repo root
            fid = fopen(cfgListPath);
            textscan(fid,'%s',6,'delimiter','\n'); % Skip header lines
            configData = textscan(fid,'%[^:]::%[^:]::%[^:\r\t\n]%*[\r\t\n]','collectoutput',true);
            configData = configData{1};
            fclose(fid);
            
            for i = 1:size(configData,1)
                if strcmp(varargin{1},configData{i,1});
                    fprintf('Initializing configuration with model %s\n',configData{i,1});
                    fprintf('"%s"\n',configData{i,3});
                    fid2 = fopen(configData{i,2});
                    execLines = textscan(fid2,'%s','delimiter','\n');
                    execLines = execLines{1};
                    commentLines = cellfun(@(x) numel(x) == 0 || (numel(x) > 0) && (x(1) == '%'),execLines,'uni',true);
                    execLines(commentLines) = [];
                    for l = 1:numel(execLines)
                        k = find(execLines{l} ~= ' ',1);
                        try
                            eval(sprintf('obj.%s;',execLines{l}(k:end)));
                        catch error
                            throw(MException('',sprintf(...
                                'mVisionConfig::mVisionConfig - Invalid custom configuration command at line %u for tag %s\n"...%s..."',...
                                l+numel(commentLines),varargin{1},execLines{l}(k:end))));
                        end
                    end
                    fclose(fid2);
                    return;
                end
            end
            throw(MException('',sprintf('mVisionConfig::mVisionConfig - Requested configuration tag "%s" is unknown.',varargin{1})));
        end % Constructor
        
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
            cleanConfig.EIRP = obj.EIRightPoints;
            cleanConfig.EISp = obj.EISpikesToAverage;
            cleanConfig.EInThreads = obj.EInThreads;
            
            cleanConfig.resamplePitch = obj.EIResamplePitch;
            cleanConfig.eiThrW = obj.EIMergeThresholdWithin;
            cleanConfig.eiThrG = obj.EIMergeThresholdGlobal;
            
            cleanConfig.globMinWin = obj.EIGlobalMinWindow;
        end
        
        function pools = buildElPools(obj,nElectrodes)
            cuts = round((0:obj.nWorkers) .* nElectrodes ./ obj.nWorkers) + 1;
            pools = cell(obj.nWorkers,1);
            for p = 1:obj.nWorkers
                pools{p} = cuts(p):(cuts(p+1)-1);
            end
        end
        
        function computeConfig = getComputeConfig(obj)
            computeConfig.eiRaw = obj.computeRawEIs;
            computeConfig.staRaw = obj.computeRawSTAs;
            computeConfig.nRaw = obj.keepNeuronsRaw;
            
            computeConfig.ei = obj.computeFinalEIs;
            computeConfig.sta = obj.computeFinalSTAs;
            computeConfig.params = obj.computeFinalParams;
        end
        
        % Generate the system call strings for .ei(-raw) computation
        % Requires EI calculation
        function commands = stringifyEICommand(obj, rawDataPath,saveFolder,datasetName)
            start = ['java -cp "./vision/',pathsep,...
                './vision/Vision.jar" edu.ucsc.neurobiology.vision.calculations.CalculationManager'];
            commands = {...
                sprintf('%s "Electrophysiological Imaging Fast" %s %s %5.3f %u %u %u %u',...
                start, saveFolder, rawDataPath, 0.1, obj.EILeftPoints, obj.EIRightPoints,...
                obj.EISpikesToAverage, obj.EInThreads) };
        end
        
        % Generate the system call strings for .sta(-raw) computation
        % STA computation requires
        % Generate globals file
        % Copy raw data header to globals
        % Make white noise movie -----------> TBD directly in the matlab-JVM workspace
        % Calculate auxiliary parameters --|
        % STA Calculation
        function commands = stringifySTACommand(obj, rawDataPath, saveFolder, datasetName)
            start = ['java -cp "./vision/',pathsep,...
                './vision/Vision.jar" edu.ucsc.neurobiology.vision.calculations.CalculationManager'];
            if isunix
                tail = '> /dev/null';
            else
                tail = '> NUL';
            end
            moviePath = [saveFolder,filesep,datasetName,'.movie'];
            commands = {...
                sprintf('%s "Generate Globals Files" %s false %s',...
                start, saveFolder, tail),...
                sprintf('%s "Copy Raw Data Header to Globals" %s %s/%s.globals %s',...
                start, rawDataPath, saveFolder, datasetName, tail),...
                sprintf('%s "STA Calculation Parallel" %s %s %u %d %u %u %u 1 0 1.0 false false',...
                start, saveFolder, moviePath, obj.STADepth, obj.STAOffset, obj.STAnSpikes,...
                obj.STACalcAtOnce, obj.STAnThreads) };
        end
        
        % Generate the system call strings for .params computation
        % Recompute a final globals and a copy raw data header to globals
        % As for unknown reasons the EI display is very erroneous and it seems
        % to fix it
        function commands = stringifyParamsCommand(obj, rawDataPath, saveFolder, datasetName)
            start = ['java -cp "./vision/',pathsep,...
                './vision/Vision.jar" edu.ucsc.neurobiology.vision.calculations.CalculationManager'];
            if isunix
                tail = '> /dev/null';
            else
                tail = '> NUL';
            end
            commands = {...
                sprintf('%s "Make Parameters File" %s %u true true 1.0 5.0 1 false true 3.0 0 0.02551 -40.0 true 100.0 0.5 false x 30 false x x false x false x true',...
                start, saveFolder, obj.paramsnThreads),...
                sprintf('%s "Generate Globals Files" %s false %s',...
                start, saveFolder, tail),...
                sprintf('%s "Copy Raw Data Header to Globals" %s %s/%s.globals %s',...
                start, rawDataPath, saveFolder, datasetName, tail)};
        end
    end % methods
end

