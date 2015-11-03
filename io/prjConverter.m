function prjConverter(dataPath, saveFolder, datasetName)
    % Converts the prj.mat file of a dataset into a vision
    % compatible .prj
    
    prjFilePath = [saveFolder,filesep,datasetName,'.prj'];
    matPrjFilePath = [saveFolder,filesep,datasetName,'.prj.mat'];
    matSpikeFilePath = [saveFolder,filesep,datasetName,'.spikes.mat'];
    matCovFilePath = [saveFolder,filesep,datasetName,'.cov.mat'];
    load(matSpikeFilePath,'nSamples');
    
    % ------------- DO NOT TRY TO UNDERSTAND THIS SECTION -----------------
    % --------- there are workarounds, but it's just all pain -------------
    % Get the bin file header
    rawDataFile = edu.ucsc.neurobiology.vision.io.RawDataFile(dataPath);
    rawDataHeader = rawDataFile.getHeader();
    rawDataFile.close();
    
    % See visionParams...
    VISION_VERSION = 8001010;
    VERSION = 32;
    
    % Instantiate a new Vision header
    visionHeader = edu.ucsc.neurobiology.vision.io.VisionHeader();
    
    % Reasonable defaults
    MIN_SPIKES = 100;
    MAX_CONTAM = 1;
    
    % Fill in the Vision header object
    visionHeader.magic = MAGIC;
    visionHeader.headerVersion = 1;
    visionHeader.version = VERSION;
    visionHeader.meanTimeConstant = -1;
    visionHeader.threshold = -1;
    visionHeader.arrayID = rawDataHeader.getArrayID();
    visionHeader.nSamples = nSamples;
    visionHeader.samplingFrequency = rawDataHeader.getSamplingFrequency();
    visionHeader.visionVersion = VISION_VERSION;
    cfg = mVisionConfig();
    cfg = cfg.getProjConfig();
    visionHeader.nDimensions = cfg.nDims;
    visionHeader.minNeuronSpikes = MIN_SPIKES;
    visionHeader.maxContamination = MAX_CONTAM;
    %----------------------------------------------------------------------
    
    load(matSpikeFilePath,'ttlTimes');
    load(matCovFilePath,'totSpikes');
    load(matPrjFilePath,'spikeTimes')
    
    prjFile = edu.ucsc.neurobiology.vision.io.ProjectionsFile(prjFilePath,...
        visionHeader, totSpikes, ttlTimes);
    for el = 2:size(totSpikes,1)
        load(matPrjFilePath,sprintf('projSpikes%u',el));
        eval(sprintf('prjFile.saveDataFromMatlab(el-1,spikeTimes{el},projSpikes%u);',el));
        eval(sprintf('projSpikes%u = [];',el));
    end
    
    prjFile.close();
end