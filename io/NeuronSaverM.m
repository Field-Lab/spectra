classdef NeuronSaverM < handle
    %NEURONSAVERM Utility class to save neurons into vision compatible format
    % Uses initVisionNeuronFile to properly create a NeuronFile
    % Wraps the java methods
    % Allows to push neurons in the newly created file
    
    properties
        nElectrodes
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
        function obj = NeuronSaverM(dataPath,saveFolder,datasetName)
            % Build file path
            obj.neuronFilePath = [saveFolder,filesep,datasetName,'.neurons-raw'];
            spikeFilePath = [saveFolder,filesep,datasetName,'.spikes.mat'];
            
            load(spikeFilePath,'ttlTimes'); % Load ttl times
            
            % Creating data source and catching nElectrodes
            dataSource = DataFileUpsampler(dataPath);
            obj.nElectrodes = dataSource.nElectrodes;
           
            % Adjusting bin path to catch header properly
            if strcmp(dataPath((end-3):end),'.bin')
                binPath = dataPath;
            else
                files = dir([dataPath,filesep,datasetName,'*.bin']);
                binPath = [dataPath,filesep,files(1).name];
            end
            
            % Initialize Neuron File
            obj.neuronFile = initVisionNeuronFile(binPath, obj.neuronFilePath, ttlTimes);
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
            obj.neuronFile.addNeuron(el-1,neurID,spikeTimes,numel(spikeTimes));
        end
        
        % Closes the Neuron File
        % Necessary
        function close(obj)
            obj.neuronFile.close();
        end
       
        % Destructor
        % Will be called by the garbage collector and will close the file if not done yet.
        function delete(obj)
            try
                obj.neuronFile.close();
            catch
            end
        end % Destructor
        
    end
    
end
