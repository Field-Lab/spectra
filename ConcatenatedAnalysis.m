function ConcatenatedAnalysis( dataCommand, saveRoot )
    %CONCATENATEDANALYSIS Summary of this function goes here
    %   Detailed explanation goes here
    
    %% Load mVision configuration - start parallel pool if neede
    config = mVisionConfig();
    parConfig = config.getParConfig();
    
    %% Parse input string
    javaaddpath ./vision/Vision.jar
    addpath ./io
    [datasets, timeCommands] = DataFileUpsampler.getDatasets(dataCommand)
    saveFoldersAndName = cell(size(datasets));
    
		% Open a parpool for parallel dataset processing
    if numel(datasets) < parConfig.nWorkers
    %    parpool(numel(datasets));
    else
    %    parpool(parConfig.nWorkers);
    end
    % Process all datasets noise + spikes + cov
    for d = 1:numel(datasets)
        [~,datasetName,~] = fileparts(datasets{d});
        saveFolderSub = [saveRoot,filesep,datasetName];
        saveFolderAndName{d} = [saveFolderSub,filesep,datasetName];
        % mVision(datasets{d}, saveFolderSub, timeCommands{d}, '', 'noisetocov','all');
    end

    % Build global spike file and cov file
    addpath ./util
    mergeCovAndSpikes(saveRoot, saveFolderAndName);
    
    % Projections for all datasets
    for d = 1:numel(datasets)
        [~,datasetName,~] = fileparts(datasets{d});
        saveFolderSub = [saveRoot,filesep,datasetName];
        % mVision(datasets{d}, saveFolderSub, timeCommands{d}, '', [0 0 0 1 0 0 0],'all');
    end
    
    mergePrj(saveRoot, saveFolderAndName);
    
    % Prj + Clusters for global dataset
    % manage parpool
		x = gcp;
		if x.NumWorkers < parConfig.nWorkers
	%			delete(gcp);
	%			parpool(parConfig.nWorkers);
		end
		mVision([filesep,'concat'], saveRoot, '', '', 'prjtoneurons','all');
    
    % Build .neurons for each file from global neurons
    % just load final neurons, and make time intersect
    % delete(gcp);
end

% ConcatenatedAnalysis('/Volumes/Archive/2005-04-26-0/data002(1780-)-data003(-20)','/Volumes/Lab/Projects/spikesorting/mvision/outputConcatTest')
