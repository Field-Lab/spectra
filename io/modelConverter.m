function modelConverter( dataPath, saveFolder, datasetName )
    % Converts the model.mat file of a dataset into a vision
    % compatible .model (hopefully)
    
    global GLOBAL_CONFIG
    covConfig = GLOBAL_CONFIG.getCovConfig();
    
    modelFilePath = [saveFolder,filesep,datasetName,'.model'];
    matModelFilePath = [saveFolder,filesep,datasetName,'.model.mat'];
    matSpikeFilePath = [saveFolder,filesep,datasetName,'.spikes.mat'];
    matPrjFilePath = [saveFolder,filesep,datasetName,'.prj.mat'];
    load(matSpikeFilePath,'nSamples');
    
    datasource = DataFileUpsampler(dataPath);
    
    [adj,~] = catchAdjWJava(datasource,covConfig.electrodeUsage);
    datasource = [];
    
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
    % visionHeader.magic = MAGIC;
    visionHeader.headerVersion = 1;
    visionHeader.version = VERSION;
    visionHeader.meanTimeConstant = -1;
    visionHeader.threshold = -1;
    visionHeader.arrayID = rawDataHeader.getArrayID();
    visionHeader.nSamples = nSamples;
    visionHeader.samplingFrequency = rawDataHeader.getSamplingFrequency();
    visionHeader.visionVersion = VISION_VERSION;
    
    projConfig = GLOBAL_CONFIG.getProjConfig();
    visionHeader.nDimensions = projConfig.nDims;
    visionHeader.minNeuronSpikes = MIN_SPIKES;
    visionHeader.maxContamination = MAX_CONTAM;
    visionHeader.nEMSpikes = 3000;
        %----------------------------------------------------------------------
    
    load(matModelFilePath);
    emptyCells = cellfun(@isempty, clusterParams);
        noEmpties = clusterParams;
        noEmpties(emptyCells) = [];
        nTotalClust = sum(cellfun(@(x) x.numClusters,noEmpties,'uni',true));
    load(matPrjFilePath,'eigenVectors');
    
    modelFile = edu.ucsc.neurobiology.vision.io.ClusteringModelFile(modelFilePath,...
        nTotalClust + 1000, visionHeader);
    
    for el = 2:numel(clusterParams)
                    mw = edu.ucsc.neurobiology.vision.io.EMModelWrapper;
            m = mw.model;
                    m.extractionID = el-1;
        if isempty(clusterParams{el})
            m.nClusters = 1;
                    m.neuronIndex = 0;
            m.neuronID = NeuronSaverM.getIDs(el,1);
            m.cleaningLevel = 1;
        
            m.electrodes = adj{el}-1;
        
            m.threshold = 0;
            m.nDimensions = 5;
            m.eigenvectors = zeros(5);
            m.nGaussians = 1;
            m.means = zeros(1,5);
            m.covariances = zeros(1,5);
            m.probability = 1;
                else
            m.nClusters = clusterParams{el}.numClusters;
                    m.neuronIndex = 0:(clusterParams{el}.numClusters - 1);
            m.neuronID = NeuronSaverM.getIDs(el,1:clusterParams{el}.numClusters);
            m.cleaningLevel = 1;
        
            m.electrodes = adj{el}-1;
        
            m.threshold = 0;
            m.nDimensions = size(clusterParams{el}.centroids,2);
            m.eigenvectors = eigenVectors{el}';
            m.nGaussians = clusterParams{el}.numClusters;
            m.means = clusterParams{el}.centroids;
            m.covariances = clusterParams{el}.covariances;
            m.probability = clusterParams{el}.mixFrac;
        end

                modelFile.addExtraction(m);
    end

        modelFile.close();

end

