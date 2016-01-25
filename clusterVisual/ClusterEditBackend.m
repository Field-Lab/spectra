classdef ClusterEditBackend < handle
    %CLUSTEREDITBACKEND Handles persistently the backend of the cluster visualizer
    
    properties (GetAccess = public, SetAccess = immutable)
        prj % struct - {bool exists ; String path}
        model % struct - {bool exists ; String path}
        neurons % struct - {bool exists ; String path}
        
        config
        displayPoints = 8000
    end
    
    properties(GetAccess = public,SetAccess = protected)
        nElectrodes
        
        neuronEls
        neuronClusters
        neuronIDs
        neuronSpikeTimes
        
        elSpikeTimes
        
        elLoaded
        prjLoaded % only handles 1 loaded prj. Caching effect TODO
        
        % Front end friendly data
        isDataReady
        
        nClusters
        displayIDs
        % contaminationValues
        % ACF
        spikeTrains % for spike rate
        prjTrains
        % cleaning status
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
        end
        
        function msg = loadEl(obj,el)
            if el == obj.elLoaded
                msg = sprintf('Electrode %u already loaded.',el);
                return
            end
            if (el <= 1) || (el > obj.nElectrodes)
                msg = sprintf('Electrode number %u invalid, nothing done.',el);
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
            
            neuronIndices = find(obj.neuronEls == el);
            reductionFrac = min(1,obj.displayPoints ./ size(obj.prjLoaded,1));
            
            for c = 1:obj.nClusters
                [~,~,indices] = intersect(obj.neuronSpikeTimes{neuronIndices(c)},obj.elSpikeTimes{obj.elLoaded});
                obj.spikeTrains{c} = obj.neuronSpikeTimes{neuronIndices(c)};
                if reductionFrac < 1
                    subind = randsample(numel(obj.spikeTrains{c}),floor(reductionFrac*numel(obj.spikeTrains{c})));
                    obj.prjTrains{c} = obj.prjLoaded(indices(subind),:);
                end
            end
            
            obj.isDataReady = true;
            msg = sprintf('Displaying electrode %u.',el);
        end
        
        function msg = loadID(obj,ID)
            rowNum = find(obj.neuronIDs == ID);
            if numel(rowNum) ~= 1
                msg = sprintf('Neuron ID invalid, nothing done.');
                return;
            else
                el = obj.neuronEls(rowNum);
                msg = obj.loadEl(el);
            end
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

