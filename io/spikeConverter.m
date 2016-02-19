function spikeConverter(dataPath, saveFolder, datasetName)
    % Converts the spikes.mat file of a dataset into a vision
    % compatible .spikes
    global GLOBAL_CONFIG
    
    spikeFilePath = [saveFolder,filesep,datasetName,'.spikes'];
    matSpikeFilePath = [saveFolder,filesep,datasetName,'.spikes.mat'];
            
    load(matSpikeFilePath);
    % This loads spikeSave, ttlTimes and nSamples
    
    cfg = GLOBAL_CONFIG.getSpikeConfig();
    threshold = cfg.spikeThreshold;
    meanTimeConstant = cfg.meanTimeConstant;
    
    datasource = DataFileUpsampler(dataPath);
    samplingRate = datasource.samplingRate;
    hdr = datasource.rawDataFile.getHeader();
    arrayID = hdr.getArrayID();
    
    spikeFile = edu.ucsc.neurobiology.vision.io.SpikeFile(spikeFilePath,...
        arrayID, meanTimeConstant, threshold, nSamples, samplingRate);
    
    spikeFile.pushMatlabSpikes(spikeSave(:,1),spikeSave(:,2)-1);
    spikeFile.close();
end

