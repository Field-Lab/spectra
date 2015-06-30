classdef NeuronSaverM < handle
    %NEURONSAVERM Utility class to save neurons into vision compatible format
    % Uses initVisionNeuronFile.m to properly create a edu.ucsc.neurobiology.vision.io.NeuronFile
    % Wraps the java method
    % Allows to push neurons
    
    properties
        nElectrodes
        neuronFilePath
        neuronFile
    end
    
    methods
        function obj = NeuronSaverM(dataPath,saveFolder,datasetName)
            obj.neuronFilePath = [saveFolder,filesep,datasetName,'.neurons'];
            spikeFilePath = [saveFolder,filesep,datasetName,'.spikes.mat'];
            
            load(spikeFilePath,'ttlTimes');
            
            %% Creating data source and catching nElectrodes
            dataSource = DataFileUpsampler(dataPath);
            
            header = dataSource.rawDataFile.getHeader();
            packedArrayID = int32(header.getArrayID());
           
            electrodeMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(packedArrayID);
            obj.nElectrodes = electrodeMap.getNumberOfElectrodes();
            
            %% Adjusting bin path to catch header properly
            if strcmp(dataPath((end-3):end),'.bin')
                binPath = dataPath;
            else
                files = dir([dataPath,filesep,datasetName,'*.bin']);
                binPath = [dataPath,filesep,files(1).name];
            end
            
            obj.neuronFile = initVisionNeuronFile(binPath, obj.neuronFilePath, ttlTimes);
        end
        
        % Input: Matlab electrode and cluster numbers
        % Pass num-1 to java
        function neurID = getNeuronID(obj, elNum, clustNum)
            
            neurID = obj.neuronFile.getNeuronID(elNum-1, clustNum-1);
        end
        
        % Input: Matlab el number
        % Pass num-1 to java
        function addNeuron(obj, el, neurID, spikeTimes)
            obj.neuronFile.addNeuron(el-1,neurID,spikeTimes,numel(spikeTimes));
        end
        
        function close(obj)
            obj.neuronFile.close();
        end
    end
    
end

