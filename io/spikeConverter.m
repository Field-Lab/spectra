function spikeConverter(dataPath, saveFolder, datasetName)
    % Converts the spikes.mat file of a dataset into a vision
    % compatible .spikes
    
    spikeFilePath = [saveFolder,filesep,datasetName,'.spikes'];
    matSpikeFilePath = [saveFolder,filesep,datasetName,'.spikes.mat'];
            
    load(matSpikeFilePath);
    % This loads spikeSave, ttlTimes and nSamples
    
    cfg = mVisionConfig();
    cfg = cfg.getSpikeConfig();
    threshold = cfg.spikeThreshold;
    meanTimeConstant = cfg.meanTimeConstant;
    
    datasource = DataFileUpsampler(dataPath);
    samplingRate = datasource.samplingRate;
    hdr = datasource.RawDataFile.getHeader();
    arrayID = hdr.getArrayID();
    
    spikeFile = edu.ucsc.neurobiology.vision.io.SpikeFile(spikeFilePath,...
        arrayID, meanTimeConstant, threshold, nSamples, samplingRate);
    
    spikeFile.pushMatlabSpikes(bsxfun(@minus,spikeSave,[0 1]));
    spikeFile.close();
end

