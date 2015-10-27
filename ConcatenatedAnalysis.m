function ConcatenatedAnalysis( dataCommand, saveRoot )
    %CONCATENATEDANALYSIS Summary of this function goes here
    %   Detailed explanation goes here
    
    %% Load mVision configuration - start parallel pool if neede
    config = mVisionConfig();
    parConfig = config.getParConfig();
    
    %% Parse input string

    [datasets, timeCommands] = DataFileUpsampler.getDatasets(dataCommand);
    % Open a parpool for parallel dataset processing
    if numel(datasets < parConfig.nWorkers);
        parpool(numel(datasets));
    else
        parpool(parConfig.nWorkers);
    end
    % Process all datasets noise + spikes + cov
    parfor d = 1:numel(datasets)
        [~,datasetName,~] = fileparts(datasets{d});
        saveFolderSub = [saveRoot,filesep,datasetName];
        mVision(datasets{d}, saveFolderSub, timeCommands{d}, '', 'noisetocov','all');
    end
    
    % Build global spike file and cov file
    addpath ./util
    mergeConcatFiles(saveRoot, datasets, timeCommands);
    
    % Prj + Clusters for global dataset
    mvision('', saveRoot, '', 'prjtoneurons','all');
    
    % Build .neurons for each file from global neurons
    % just load final neurons, and make time intersect
    delete(gcp);
end

