classdef ClusterEditBackend < handle
    %CLUSTEREDITBACKEND Handles persistently the backend of the cluster visualizer
    
    properties (GetAccess = public, SetAccess = immutable)
        analysisPath
        
        prj % struct - {bool exists ; String path}
        model % struct - {bool exists ; String path}
        neurons % struct - {bool exists ; String path}
        
        config
        
        % Java related stuff for ei and sta display
        eiFile
        staFile
        globalsFile
        
        arrayID
        electrodeMap
    end
    
    properties(GetAccess = public,SetAccess = protected)
        nNeurons
        nElectrodes
        nSamples
        
        neuronEls
        neuronClusters
        neuronIDs
        neuronStatuses
        classification
        
        elLoaded
        prjLoaded % only handles 1 loaded prj. Caching effect TODO
        
        % Front end friendly data
        isDataReady
        
        nClusters
        displayIDs
        contaminationValues
        spikeCounts
        statusRaw
        comment
        % ACF
        spikeTrains % for spike rate
        spikeTrainCorr
        prjTrains
        eisLoaded
        EIdistMatrix
        stasLoaded
        % cleaning status
    end
    
    properties % public setting access
        statusBarHandle
    end
    
    methods
        % Constructor
        % Builds backend object defined in this file for the cluster visualizer
        % Input
        %       analysisPath: path to analysis folder
        %       OPTIONAL (none or all three)
        %       pathToPrj: full path of .prj.mat file to use
        %       pathToNeurRaw: full path to .neurons.mat to use
        %       pathToModel: full path to .model.mat to use (unused so far)
        %
        % One or more of the optional arguments can be set to empty string ''
        % in which case constructor looks for a file with correct extension
        % in the analysis path.
        function obj = ClusterEditBackend(analysisPath,varargin)
            obj.analysisPath = analysisPath;
            obj.config = mVisionConfig();
            
            % Argument check
            validateattributes(analysisPath,{'char'},{},'','Analysis Path',1);
            if ~ (exist(analysisPath,'file') == 7)
                throw(MException('','ClusterEditBackend:ClusterEditBackend - Analysis folder does not exist'));
            end
            
            if nargin > 1
                narginchk(4,4);
                validateattributes(varargin{1},{'char'},{},'','projections file',2);
                validateattributes(varargin{2},{'char'},{},'','neuron file',3);
                validateattributes(varargin{3},{'char'},{},'','model file',4);
            end
            
            % nSamples - do that clean
            files = dir([analysisPath,filesep,'*.spikes.mat']);
            spikePath = [analysisPath,filesep,files(1).name];
            load(spikePath,'nSamples');
            obj.nSamples = nSamples;
            
            % Look for a .prj.mat, a .neurons.mat and a .model.mat
            % Projections file
            if nargin == 1 || (nargin == 4 && strcmp(varargin{1},'')) % No forced file as optional argument
                files = dir([analysisPath,filesep,'*.prj.mat']);
            else % File manually selected as optional argument
                files = dir(varargin{1});
            end
            if numel(files) == 1
                obj.prj.exist = true;
                obj.prj.path = [analysisPath,filesep,files(1).name];
            else if numel(files) == 0
                    throw(MException('','ClusterEditBackend:ClusterEditBackend - No projections file found'));
                else
                    throw(MException('','ClusterEditBackend:ClusterEditBackend - Multiple projections files found\nPlease use extended constructor call\nneuronViewer(analysisPath,projectionsFile,neuronsRawFile,modelFile)'));
                end
            end
            
            % Neurons file
            if nargin == 1 || (nargin == 4 && numel(varargin{2}) == 0)
                files = dir([analysisPath,filesep,'*.neurons.mat']);
            else
                files = dir(varargin{2});
            end
            if numel(files) == 1
                obj.neurons.exist = true;
                obj.neurons.path = [analysisPath,filesep,files(1).name];
            else if numel(files) == 0
                    throw(MException('','ClusterEditBackend:ClusterEditBackend - No neurons-raw file found'));
                else
                    throw(MException('','ClusterEditBackend:ClusterEditBackend - Multiple neurons-raw files found\nPlease use extended constructor call\nneuronViewer(analysisPath,projectionsFile,neuronsRawFile,modelFile)'));
                end
            end
            
            % (Optional) model file
            if nargin == 1 || (nargin == 4 && numel(varargin{3}) == 0)
                files = dir([analysisPath,filesep,'*.model.mat']);
            else
                files = dir(varargin{3});
            end
            if numel(files) == 1
                obj.model.exist = true;
                obj.model.path = [analysisPath,filesep,files(1).name];
            else if numel(files) == 0
                    fprintf('ClusterEditBackend:ClusterEditBackend - No model file found, skipping\n');
                    obj.model.exist = false;
                else
                    throw(MException('','ClusterEditBackend:ClusterEditBackend - Multiple model files found\nPlease use extended constructor call\nneuronViewer(analysisPath,projectionsFile,neuronsRawFile,modelFile)'));
                end
            end
            % All files now referenced and checked
            % Partial loading for neuronSpikeTimes and elSpikeTimes
            obj.prj.matfile = matfile(obj.prj.path);
            obj.neurons.matfile = matfile(obj.neurons.path);
            
            % Load everything useful and RAM OK
            % Generate metadata
            
            load(obj.neurons.path,'neuronEls','neuronClusters'); % Loads neuronClusters, neuronEls
            obj.neuronEls = neuronEls;
            obj.neuronClusters = neuronClusters;
            obj.neuronIDs = ClusterEditBackend.getIDs(obj.neuronEls, obj.neuronClusters);
            
            obj.nElectrodes = size(obj.prj.matfile.spikeTimes,1);
            obj.nNeurons = size(obj.neuronEls,1);
            
            % CleanPattern - do that clean
            cleanPatternPath = [analysisPath,filesep,'cleanPattern.mat'];
            obj.neuronStatuses = zeros(obj.nNeurons,2);
            if exist(cleanPatternPath)
                load(cleanPatternPath);
                IDsDuplicatesRemoved = ClusterEditBackend.shortenDuplicatesPath(IDsDuplicatesRemoved);
                [~,i,~] = intersect(obj.neuronIDs,IDsRemovedAtContam);
                obj.neuronStatuses(i,1) = 1; % 1 - removed at contam
                [~,i,j] = intersect(obj.neuronIDs,IDsMerged(:,2));
                obj.neuronStatuses(i,1) = 2; % 2 - merged
                obj.neuronStatuses(i,2) = IDsMerged(j,1);
                if size(IDsDuplicatesRemoved,2) == 1 % old format, no track
                    [~,i,~] = intersect(obj.neuronIDs,IDsDuplicatesRemoved(:,1));
                    obj.neuronStatuses(i,1) = 3; % 3 - duplicate
                    obj.neuronStatuses(i,2) = -1;
                else
                    [~,i,j] = intersect(obj.neuronIDs,IDsDuplicatesRemoved(:,2));
                    obj.neuronStatuses(i,1) = 3; % 3 - duplicate
                    obj.neuronStatuses(i,2) = IDsDuplicatesRemoved(j,1);
                end
            end
            
            % Initialize EI
            files = dir([analysisPath,filesep,'*.ei']);
            if numel(files) > 0
                eiPath = [analysisPath,filesep,files(1).name];
                obj.eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(eiPath);
            else
                obj.eiFile = []; % Tag for existence of ei File
                fprintf('ClusterEditBackend:ClusterEditBackend - No EI file found, skipping\n');
            end
            
            % Initialize STA
            files = dir([analysisPath,filesep,'*.sta']);
            if numel(files) > 0
                staPath = [analysisPath,filesep,files(1).name];
                obj.staFile = edu.ucsc.neurobiology.vision.io.STAFile(staPath);
            else
                obj.staFile = []; % Tag for existence of ei File
                fprintf('ClusterEditBackend:ClusterEditBackend - No STA file found, skipping\n');
            end
            
            % Grab the globals file
            files = dir([analysisPath,filesep,'*.globals']);
            if numel(files) > 0
                globalsPath = [analysisPath,filesep,files(1).name];
                obj.globalsFile = edu.ucsc.neurobiology.vision.io.chunk.GlobalsFile(globalsPath);
            else
                obj.globalsFile = [];
                obj.staFile = []; % Tag for existence of ei File
                fprintf('ClusterEditBackend:ClusterEditBackend - No Globals file found, STAs Unavailable\n');
            end
            
            % Try to find a classification file
            obj.classification = cell(obj.nNeurons,1);
            files = dir([analysisPath,filesep,'*.txt']);
            if numel(files) == 1
                fid = fopen([analysisPath,filesep,files(1).name]);
                classesRaw = textscan(fid, '%u All/%s', 'delimiter', '\n');
                [~,pos,pos2] = intersect(obj.neuronIDs,classesRaw{1});
                obj.classification(pos) = classesRaw{2}(pos2);
                fclose(fid);
            else
                fprintf('ClusterEditBackend:ClusterEditBackend - Can''t find classification .txt file\n.');
            end
            
            % electrode map
            if numel(obj.eiFile) > 0
                obj.arrayID = obj.eiFile.arrayID;
                obj.electrodeMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(obj.arrayID);
            end
        end
        
        function returnStatus = loadEl(obj,el)
            if el == obj.elLoaded
                obj.statusBarHandle.String = sprintf('Electrode %u already loaded.',el-1);
                returnStatus = 1;
                return;
            end
            if (el <= 1) || (el > obj.nElectrodes)
                obj.statusBarHandle.String = sprintf('Electrode number %u invalid, nothing done.',el-1);
                returnStatus = 1;
                return;
            else
                obj.elLoaded = el;
                eval(sprintf('load(''%s'',''projSpikes%u'');',obj.prj.path,el));
                eval(sprintf('obj.prjLoaded = projSpikes%u;',el));
            end
            
            neuronIndices = find(obj.neuronEls == el);
            
            obj.displayIDs = obj.neuronIDs(neuronIndices);
            obj.nClusters = numel(obj.displayIDs);
            obj.spikeTrains = cell(obj.nClusters,1);
            obj.prjTrains = cell(obj.nClusters,1);
            obj.contaminationValues = zeros(obj.nClusters,1);
            obj.spikeCounts = zeros(obj.nClusters,1);
            obj.statusRaw = obj.neuronStatuses(neuronIndices,:);
            obj.comment = obj.classification(neuronIndices);
            
            elSpikeTimes = obj.prj.matfile.spikeTimes(obj.elLoaded,1);
            elSpikeTimes = elSpikeTimes{1};
            
            if numel(neuronIndices) > 0
                obj.spikeTrains = obj.neurons.matfile.neuronSpikeTimes(neuronIndices,1);
            else
                obj.spikeTrains = cell(0,1);
            end
            obj.spikeCounts = cellfun(@numel, obj.spikeTrains,'uni',true);
            
            for c = 1:obj.nClusters
                [~,~,indices] = intersect(...
                    obj.spikeTrains{c},...
                    elSpikeTimes);
                obj.prjTrains{c} = obj.prjLoaded(indices,:);
                
                obj.contaminationValues(c) = ...
                    edu.ucsc.neurobiology.vision.anf.NeuronCleaning.getContam(obj.spikeTrains{c},int32(obj.nSamples));
            end
            
            if numel(obj.eiFile) > 0
                obj.loadEI();
                obj.loadEIDistances();
            end
            if numel(obj.staFile) > 0
                obj.loadSTA();
            end
            obj.loadCorrelations();
            
            obj.isDataReady = true;
            returnStatus = 0;
        end
        
        function loadEI(obj)
            obj.eisLoaded = cell(obj.nClusters,1);
            for c = 1:obj.nClusters
                try
                    obj.eisLoaded{c} = obj.eiFile.getImage(obj.displayIDs(c));
                catch
                    obj.eisLoaded{c} = [];
                end     
            end
        end
        
        function loadEIDistances(obj)
            % Purpose of this function is to select relevant electrodes as in merging,
            % Realign EIs and compute the matrix of distances.
            % The process is meant to be identical as merging analysis in duplicate removal
            cfg = mVisionConfig();
            cleanConfig = cfg.getCleanConfig();
            
            eiSize = cleanConfig.EILP + cleanConfig.EIRP + 1;
            minWindowSize = cleanConfig.globMinWin(2) - cleanConfig.globMinWin(1);
            samplingVal = cleanConfig.globMinWin(1):cleanConfig.resamplePitch:cleanConfig.globMinWin(2);
            basicValues = 1:(eiSize - minWindowSize);
            
            [adjacent2,~] = catchAdjWJava( obj.eiFile, 2);
            
            trace = zeros(obj.nClusters, (eiSize - minWindowSize) * numel(adjacent2{obj.elLoaded}));
            
            % EI upsampling handling
            traceTemp = zeros(eiSize - minWindowSize, numel(adjacent2{obj.elLoaded}));
            hasEI = cellfun(@(x) numel(x) > 0,obj.eisLoaded);
            for n = 1:obj.nClusters
                if ~hasEI(n)
                    continue
                end
                % Realignment procedure
                for neighborElInd = 1:numel(adjacent2{obj.elLoaded});
                    neighborEl = adjacent2{obj.elLoaded}(neighborElInd);
                    
                    interpolant = griddedInterpolant(1:eiSize,squeeze(obj.eisLoaded{n}(1,neighborEl,:)),'spline');
                    [~, offset] = min(interpolant(samplingVal));
                    traceTemp(:,neighborElInd) = interpolant(basicValues + (offset - 1) * cleanConfig.resamplePitch);
                end
                trace(n,:) = traceTemp(:)';
            end % n
            % L2 normalization
            trace(hasEI,:) = bsxfun(@rdivide,trace(hasEI,:),sqrt(sum(trace(hasEI,:).^2,2)));
            
            % Distance computation
            [~,dists] = returnL2MergeClasses(trace(hasEI,:), 0);
            obj.EIdistMatrix = nan(obj.nClusters);
            obj.EIdistMatrix(hasEI,hasEI) = dists;
        end
        
        function loadCorrelations(obj)
            visionCoincidenceTime = 10;
            obj.spikeTrainCorr = ones(obj.nClusters);
            for i = 1:obj.nClusters
                for j = (i+1):obj.nClusters
                    obj.spikeTrainCorr(i,j) = ...
                        edu.ucsc.neurobiology.vision.anf.NeuronCleaning.getCorrVal(...
                        obj.spikeTrains{i},obj.spikeTrains{j},...
                        visionCoincidenceTime);
                    obj.spikeTrainCorr(j,i) = obj.spikeTrainCorr(i,j);
                end
            end
        end
        
        function loadSTA(obj)
            obj.stasLoaded = cell(obj.nClusters,1);
            for c = 1:obj.nClusters
                try
                    obj.stasLoaded{c} = obj.staFile.getSTA(obj.displayIDs(c));
                catch
                    obj.stasLoaded{c} = [];
                end
            end
        end
        
        function el = checkID(obj,ID)
            rowNum = find(obj.neuronIDs == ID);
            if numel(rowNum) ~= 1
                obj.statusBarHandle.String = sprintf('Neuron ID invalid, nothing done.');
                el = -1;
                return;
            end
            [el,~] = obj.getElClust(ID);
        end
    end
    
    
    
    methods(Static)
        % Statically computes neuron IDs for electrode-clusterID pairs
        % Allows to store absolute neuron information during neuron cleaning
        % Rather than the pattern of successive removals
        % See duplicateRemoval.m, variables saved in cleanPattern.mat
        function IDs = getIDs(el, clust)
            cfg = mVisionConfig();
            cfg = cfg.getSpectralConfig();
            maxClust = cfg.maxEV;
            IDs = (el-2)*maxClust + clust;
        end
        
        function [el, clust] = getElClust(ID)
            cfg = mVisionConfig();
            cfg = cfg.getSpectralConfig();
            maxClust = cfg.maxEV;
            el = floor((ID-1)./maxClust) + 2;
            clust = mod(ID-1,maxClust) + 1;
        end
        
        function IDPairs = shortenDuplicatesPath(IDPairs)
            % 1st column - ID for which it is deleted
            % 2nd column - ID deleted
            while true
                [~,i,j] = intersect(IDPairs(:,1),IDPairs(:,2));
                if numel(i) == 0
                    break
                end
                IDPairs(i,1) = IDPairs(j,1);
            end
        end
        
    end % Static methods
end % Class
