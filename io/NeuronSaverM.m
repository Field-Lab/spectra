classdef NeuronSaverM < handle
    %NEURONSAVERM Utility class to save neurons into vision compatible format
    % Uses initVisionNeuronFile to properly create a NeuronFile
    % Wraps the NeuronFile underlying java methods
    % Allows to push neurons in the newly created file
    %
    % Author -- Vincent Deo -- Stanford University -- August 27, 2015
    
    properties
        initialOffset
        neuronFilePath
        neuronFile % edu.ucsc.neurobiology.vision.io.NeuronFile object
    end % properties
    
    methods
        % Constructor
        % Initializes the wrapper NeuronSaverM object and the underlying NeuronFile
        %
        % Inputs:
        %   datapath: path to the raw data folder
        %   saveFolder: path to the analysis folder
        %   datasetName: 'data0xx' string, to reconstruct file name
        function obj = NeuronSaverM(dataPath,saveFolder,datasetName,fileExt,initialOffset)
            % Build file path
            obj.neuronFilePath = [saveFolder,filesep,datasetName,'.neurons',fileExt];
            spikeFilePath = [saveFolder,filesep,datasetName,'.spikes.mat'];
            
            load(spikeFilePath,'ttlTimes'); % Load ttl times
            
            % required to substract the initial spike time offset if the processing did not start at
            % beginning of .bin file.
            % This is because mVision uses a different convention than vision for spike labelling in case of
            % timetag starting > 0. mVision: absolute labeling, vision: 0 at start of timeTag.
            obj.initialOffset = initialOffset; % in samples
            
            % Initialize Neuron File
            obj.neuronFile = initVisionNeuronFile(dataPath, obj.neuronFilePath, ttlTimes);
        end
        
        % Gets a Neuron ID according to Vision's standards
        % Substracts 1 to el and clust arguments and passes to equivalent java method
        %
        % Input:
        %   elNum: electrode number, in Matlab format
        %   clustNum: cluster number, in Matlab format
        function neurID = getNeuronID(obj, elNum, clustNum)
            neurID = obj.neuronFile.getNeuronID(elNum-1, clustNum-1);
        end
        
        
        % Pushes a new neuron in the NeuronFile
        % by calling the equivalent java method
        %
        % Inputs:
        %   el: electrode number in Matlab format
        %   neurID: neuron ID, vision format
        %   spikeTimes: array of spike times for this neuron
        function addNeuron(obj, el, neurID, spikeTimes)
            obj.neuronFile.addNeuron(el-1,neurID,spikeTimes-obj.initialOffset,numel(spikeTimes));
        end
        
        
        % Pushes a list of neurons all at once
        function pushAllNeurons(obj, neuronEls, neuronClusters, neuronSpikeTimes)
            nNeurons = size(neuronEls,1); 
            validateattributes(neuronEls,{'numeric'},{'size',[nNeurons, 1]},'','neuronEls');
            validateattributes(neuronClusters,{'numeric'},{'size',[nNeurons, 1]},'','neuronClusters');
            validateattributes(neuronSpikeTimes,{'cell'},{'size',[nNeurons, 1]},'','neuronSpikeTimes');
            
            for n = 1:nNeurons
                el = neuronEls(n);
                nID = obj.getNeuronID(el, neuronClusters(n));
                obj.addNeuron(el, nID, neuronSpikeTimes{n});
            end
        end
            
            
        % Closes the Neuron File
        % Necessary
        function close(obj)
            obj.neuronFile.close();
        end
       
        % Destructor
        % Will be called (?) by the garbage collector and will close the file if not done yet.
        function delete(obj)
            try
                obj.neuronFile.close();
            catch
                % If already closed, catch and do nothing.
            end
        end % Destructor
        
    end
    
    methods(Static)
        % Statically computes neuron IDs for electrode-clusterID pairs
        % Allows to store absolute neuron information during neuron cleaning
        % Rather than the pattern of successive removals
        % See duplicateRemoval.m, variables saved in .clean.mat file        
        function IDs = getIDs(el, clust)
            global GLOBAL_CONFIG
            if numel(GLOBAL_CONFIG) == 0 % Default
                GLOBAL_CONFIG = mVisionConfig();
            end
            cfg = GLOBAL_CONFIG.getSpectralConfig();
            maxClust = cfg.maxEV;
            IDs = (el-2)*maxClust + clust;
        end
        
        function [el, clust] = getElClust(ID)
            global GLOBAL_CONFIG
            if numel(GLOBAL_CONFIG) == 0 % Default
                GLOBAL_CONFIG = mVisionConfig();
            end
            cfg = GLOBAL_CONFIG.getSpectralConfig();
            maxClust = cfg.maxEV;
            el = floor((ID-1)./maxClust) + 2;
            clust = mod(ID-1,maxClust) + 1;
        end
    end
end
