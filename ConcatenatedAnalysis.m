function ConcatenatedAnalysis( dataCommand, saveRoot, varargin )
    %CONCATENATEDANALYSIS
    
    %% Load mVision configuration - start parallel pool if needed
    global GLOBAL_CONFIG
    if numel(varargin) == 1
        GLOBAL_CONFIG = mVisionConfig(varargin{1});
    else % == 0
        GLOBAL_CONFIG = mVisionConfig();
    end
    parConfig = GLOBAL_CONFIG.getParConfig();
    
    %% Parse input string
    javaaddpath ./vision/Vision.jar
    addpath(genpath(['.',filesep]));
    
    [datasets, timeCommands] = DataFileUpsampler.getDatasets(dataCommand)
    saveFolderSub = cell(size(datasets)); %stores saveRoot/dataxxx
    saveFoldersAndName = cell(size(datasets)); % stores saveRoot/dataxxx/dataxxx
    datasetName = cell(size(datasets)); %stores just dataxxx
    
    % Open a parpool for parallel dataset processing
    if numel(datasets) < parConfig.nWorkers
        parpool(numel(datasets));
    else
        parpool(parConfig.nWorkers);
    end
    % Process all datasets noise + spikes + cov
    % parfor
    parfor d = 1:numel(datasets)
        if isempty(timeCommands{d})
            timeCommands{d} = '';
        end
        [~,datasetName{d},~] = fileparts(datasets{d});
        saveFolderSub{d} = [saveRoot,filesep,datasetName{d}];
        saveFolderAndName{d} = [saveFolderSub{d},filesep,datasetName{d}];
        mVision(datasets{d}, saveFolderSub{d}, timeCommands{d}, '', 'noisetocov','all',varargin);
    end
    
    % Build global spike file and cov file
    addpath ./util
    mergeCovAndSpikes(saveRoot, saveFolderAndName, timeCommands);
    
    % Projections for all datasets
    % parfor
    for d = 1:numel(datasets)
        mVision(datasets{d}, saveFolderSub{d}, timeCommands{d}, '', [0 0 0 1 0 0],'all',varargin);
    end
    
    mergePrj(saveRoot, saveFolderAndName, timeCommands);
    
    % Clustering for global dataset
    % manage parpool
    x = gcp;
    if numel(x) == 0 || x.NumWorkers < parConfig.nWorkers
        delete(gcp);
        parpool(parConfig.nWorkers);
    end
    mVision([filesep,'concat'], saveRoot, '', '', [0 0 0 0 1 0],'all');
    
    concatenatedDuplicateRemoval(datasets, saveRoot, saveFolderAndName, timeCommands);
    
    modelConverter(datasets{1}, saveRoot, 'concat');
    spikeConverter(datasets{1}, saveRoot, 'concat');
    prjConverter(datasets{1}, saveRoot, 'concat');
    
    for d = 1:numel(datasets)
        % Broadcasting model file
        copyfile([saveRoot,filesep,'concat.model.mat'],[saveFolderAndName{d},'.model.mat']);
        % Conversion of files
        modelConverter(datasets{d}, saveFolderSub{d}, datasetName{d})
        spikeConverter(datasets{d}, saveFolderSub{d}, datasetName{d})
        prjConverter(datasets{d}, saveFolderSub{d}, datasetName{d})
    end
    
    delete(gcp);
end

