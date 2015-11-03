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
    saveFolderSub = cell(size(datasets)); %stores saveRoot/dataxxx
    saveFoldersAndName = cell(size(datasets)); % stores saveRoot/dataxxx/dataxxx
    
    % Open a parpool for parallel dataset processing
    if numel(datasets) < parConfig.nWorkers
        parpool(numel(datasets));
    else
        parpool(parConfig.nWorkers);
    end
    % Process all datasets noise + spikes + cov
    parfor d = 1:numel(datasets)
        [~,datasetName,~] = fileparts(datasets{d});
        saveFolderSub{d} = [saveRoot,filesep,datasetName];
        saveFolderAndName{d} = [saveFolderSub{d},filesep,datasetName];
        %mVision(datasets{d}, saveFolderSub{d}, timeCommands{d}, '', 'noisetocov','all');
    end
    
    % Build global spike file and cov file
    addpath ./util
    mergeCovAndSpikes(saveRoot, saveFolderAndName, timeCommands);
    
    % Projections for all datasets
    parfor d = 1:numel(datasets)
        %mVision(datasets{d}, saveFolderSub{d}, timeCommands{d}, '', [0 0 0 1 0 0 0],'all');
    end
    
    % mergePrj(saveRoot, saveFolderAndName, timeCommands);
    
    % Clustering for global dataset
    % manage parpool
    x = gcp;
    if x.NumWorkers < parConfig.nWorkers
        delete(gcp);
        parpool(parConfig.nWorkers);
    end
    mVision(datasets{1}, saveRoot, '', '', [0 0 0 0 1 0 0],'all');
    
    % Splitting global neuron file for each dataset
    splitNeurons(saveRoot, saveFolderAndName, timeCommands);
    
    % Apply neuron cleaning on first sub-dataset
    % Not ideal, but cannot compute concatened EIs effieciently and or relevantly
    mVision(datasets{1}, saveFolderSub{1}, timeCommands{1}, '', [0 0 0 0 0 1 0],'all');
    
    % Load cleaning pattern
    load([datasets{1},filesep,'cleanPattern.mat']);
    load([saveRoot,filesep,'concat.neurons.mat']);
    applyCleaningPattern( IDsRemovedAtContam, IDsMerged, IDsDuplicatesRemoved,...
        neuronEls, neuronClusters, neuronSpikeTimes );
    save([saveRoot,'concat.neurons.mat'],...
        'neuronEls','neuronClusters','neuronSpikeTimes');
    
    for d = 2:numel(datasets)
        load([saveFolderAndName{d},'.neurons.mat']);
        % Apply cleaning pattern
        applyCleaningPattern( IDsRemovedAtContam, IDsMerged, IDsDuplicatesRemoved,...
            neuronEls, neuronClusters, neuronSpikeTimes );
        save([saveFolderAndName{d},'.neurons.mat'],...
            'neuronEls','neuronClusters','neuronSpikeTimes'); 
    end
    
    % Now convert all relevant files to vision format...
    
% parfor mVision(datasets{d}, saveFolderSub{d}, timeCommands{d}, '', [0 0 0 0 0 0 1],'all');

    % Vision conversion of all required files
    delete(gcp);
end

