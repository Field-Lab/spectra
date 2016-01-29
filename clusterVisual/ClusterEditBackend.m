classdef ClusterEditBackend < handle
    %CLUSTEREDITBACKEND Handles persistently the backend of the cluster visualizer
    
    properties (GetAccess = public, SetAccess = immutable)
        prj % struct - {bool exists ; String path}
        model % struct - {bool exists ; String path}
        neurons % struct - {bool exists ; String path}
        
        config
    end
    
    properties(GetAccess = public,SetAccess = protected)
        nNeurons
        nElectrodes
        nSamples
        
        neuronEls
        neuronClusters
        neuronIDs
        neuronSpikeTimes
        neuronStatuses
        
        elSpikeTimes
        
        elLoaded
        prjLoaded % only handles 1 loaded prj. Caching effect TODO
        
        % Front end friendly data
        isDataReady
        
        nClusters
        displayIDs
        contaminationValues
        spikeCounts
        status
        % ACF
        spikeTrains % for spike rate
        prjTrains
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
                    fprintf('ClusterEditBackend:ClusterEditBackend - No model file found, skipping');
                    obj.model.exist = false;
                else
                    throw(MException('','ClusterEditBackend:ClusterEditBackend - Multiple model files found\nPlease use extended constructor call\nneuronViewer(analysisPath,projectionsFile,neuronsRawFile,modelFile)'));
                end
            end
            % All files now referenced and checked
            
            % Load everything useful and RAM OK
            % Generate metadata
            
            load(obj.neurons.path); % Loads neuronClusters, neuronEls, neuronSpikeTimes
            obj.neuronEls = neuronEls;
            obj.neuronClusters = neuronClusters;
            obj.neuronSpikeTimes = neuronSpikeTimes;
            obj.neuronIDs = ClusterEditBackend.getIDs(obj.neuronEls, obj.neuronClusters);
            
            load(obj.prj.path,'spikeTimes');
            obj.elSpikeTimes = spikeTimes;
            
            obj.nElectrodes = size(spikeTimes,1);
            obj.nNeurons = size(obj.neuronEls,1);
            
            % CleanPattern - do that clean
            cleanPatternPath = [analysisPath,filesep,'cleanPattern.mat'];
            load(cleanPatternPath);
            obj.neuronStatuses = zeros(obj.nNeurons,2);
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
        
        function returnStatus = loadEl(obj,el)
            if el == obj.elLoaded
                obj.statusBarHandle.String = sprintf('Electrode %u already loaded.',el);
                returnStatus = 1;
                return;
            end
            if (el <= 1) || (el > obj.nElectrodes)
                obj.statusBarHandle.String = sprintf('Electrode number %u invalid, nothing done.',el);
                returnStatus = 1;
                return;
            else
                obj.elLoaded = el;
                eval(sprintf('load(''%s'',''projSpikes%u'');',obj.prj.path,el));
                eval(sprintf('obj.prjLoaded = projSpikes%u;',el));
            end
            
            obj.displayIDs = obj.neuronIDs(obj.neuronEls == el);
            obj.nClusters = numel(obj.displayIDs);
            obj.spikeTrains = cell(obj.nClusters,1);
            obj.prjTrains = cell(obj.nClusters,1);
            obj.contaminationValues = zeros(obj.nClusters,1);
            obj.spikeCounts = zeros(obj.nClusters,1);
            obj.status = cell(obj.nClusters,1);
            
            neuronIndices = find(obj.neuronEls == el);
            
            for c = 1:obj.nClusters
                [~,~,indices] = intersect(...
                    obj.neuronSpikeTimes{neuronIndices(c)},...
                    obj.elSpikeTimes{obj.elLoaded});
                obj.spikeTrains{c} = obj.neuronSpikeTimes{neuronIndices(c)};
                obj.prjTrains{c} = obj.prjLoaded(indices,:);
                obj.spikeCounts(c) = numel(obj.spikeTrains{c});
                
                obj.contaminationValues(c) = ...
                    edu.ucsc.neurobiology.vision.anf.NeuronCleaning.getContam(obj.spikeTrains{c},int32(obj.nSamples));
                % Neuron statuses
                switch obj.neuronStatuses(neuronIndices(c),1)
                    case 0
                        obj.status{c} = 'Keep';
                    case 1
                        obj.status{c} = 'Contam / Low count';
                    case 2
                        obj.status{c} = sprintf('Merge with %u',obj.neuronStatuses(neuronIndices(c),2));
                    case 3
                        [e,~] = obj.getElClust(obj.neuronStatuses(neuronIndices(c),2));
                        obj.status{c} = sprintf('Duplicate of %i on el %u',obj.neuronStatuses(neuronIndices(c),2),e);
                end
            end
            
            obj.isDataReady = true;
            obj.statusBarHandle.String = sprintf('Displaying electrode %u.',el);
            returnStatus = 0;
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
    end
end

