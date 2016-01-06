function mVision(dataPath, saveFolder, timeCommand, tryToDo, force)
    %MVISION Main function for single dataset white noise analysis
    %
    % This script manages the successive steps of serial neuron finding
    % Input:
    %   dataSetName: path to the dataset folder
    %       ex: something ending in ...'yyyy-mm-dd-n/dataxxx'
    
    %   saveFolder: the output folder of the analysis.
    %       It is MOST recommended (ie do it) the deepest level subfolder has the same
    %       name 'dataxxx' as the dataset folder. Vision infers file names from
    %       destination folder whereas mVision does it from the source folder.
    %
    %   timeCommand: vision style time command for selection within a dataset
    %       Use empty string '' for full dataset
    %       ex: '', '(10-)', '(500-1000)','(-100)'
    %       Warning: mVision.m DOES NOT support cocatenated dataset syntaxes
    %       such as data000(1000-)-data002(-100).
    %       Use concatenatedAnalysis.m (WIP) for that purpose
    %   
    %   tryToDo: a 1x6 binary array for what subcalculations to do (see list below)
    %       0 - do not try the calculation
    %       1 - try the calculation
    %       String argument 'all' is allowed as shorthand for [1 1 1 1 1 1]
    %       If a calculation is required but the previous was not and the input
    %       files are missing, mVision will crash.
    %
    %   force: a 1x6 binary array for overwrite authorization of the calculations
    %       0 - do not overwrite, skip the calculation if its output files are found
    %       1 - overwrite, do the calculation in all cases
    %       String arguments 'all' and 'none' are allowed as shorthand for
    %       respectively all ones and all zeros.
    %
    % This function successively calls if requested or necessary the calculations:
    %       1/ Noise evaluation
    %       2/ Spike Finding
    %       3/ Covariance Calculation
    %       4/ Projections Calculations
    %       5/ Clustering
    %       6/ Neuron Cleaning and saving
    %       
    % All intermediate files are stored in convenient .mat files.
    % The final .neurons.mat file is converted to a vision compatible .neurons
    %
    % VISION POST PROCESSING
    % To run further analysis with vision, first call the calculations
    % "Generate globals file" then "Copy Raw Data Header To Globals"
    % You can the compute EIs, STAs and the parameter file.
    %
    % Author -- Vincent Deo -- Stanford University -- January 4, 2016
    
    %% SETUP
    % Add subfolders to the matlab path
    addpath(genpath(['.',filesep]));
    
    % generate repository path
    repoPath = pwd;
    
    % Java setup - using vision jar available in the repository
    visionPath = [repoPath,filesep,'vision',filesep,'Vision.jar'];
    if ~exist('edu/ucsc/neurobiology/vision/io/RawDataFile','class')
        javaaddpath(visionPath,'end')
    end
    javaaddpath(['.',filesep,'vision']);
    javaaddpath(['.',filesep,'clusterUtil']);
    javaaddpath(['.',filesep,'duplicateRemoval',filesep,'java_EI_comparison',filesep]);
    
    % USER INPUT - Set up data and output folders
    if dataPath(end) == filesep
        dataPath = dataPath(1:(end-1));
    end
    if saveFolder(end) == filesep
        saveFolder = saveFolder(1:(end-1));
    end
    if ~strcmp(dataPath,[filesep,'concat']) && ...
        ~(exist(dataPath,'file') == 2 || exist(dataPath,'file') == 7)
        % [filesep, 'concat'] is a special token here for concatenated analysis
        throw(MException('','demoScript: data source folder|file does not exist'));
    end
    
    mkdir(saveFolder);
    [~,datasetName,~] = fileparts(dataPath); % Catching dataset name as last part of saveFolder
    
    if exist(saveFolder,'file') == 2
        [~,datasetName,~] = fileparts(saveFolder);
    end
    
    % USER input - tryToDo -- won't do any task unless stated here
    nSteps = 6;
    % Steps in order:
    % --------- noise - spike - cov - prj - clust - save ----------------------
    if isa(tryToDo,'char')
        tryToDo = lower(tryToDo);
        if strcmp(tryToDo,'all')
            tryToDo = ones(1,nSteps);
        end
        if isa(tryToDo,'char')
            throw(MException('','Invalid tryToDo argument'));
        end
    else
        validateattributes(tryToDo,{'numeric'},{'row','ncols',nSteps,'binary'},'','tryToDo',4);
    end
    
    if isa(force,'char')
        force = lower(force);
        if strcmp(force,'all')
            force = ones(1,nSteps);
        end
        if strcmp(force,'none')
            force = zeros(1,nSteps);
        end
        if isa(force,'char')
            throw(MException('','Invalid force argument'));
        end
    else
        validateattributes(force,{'numeric'},{'row','ncols',nSteps,'binary'},'','force',4);
    end
    
    
    % DEBUG - additional saved file dataset name extension
    nameExt = '';
    
    totalTime = tic;
    
    %% Process noise and make a .noise file
    if tryToDo(1) &&...
            (force(1) || ~(exist([saveFolder,filesep,datasetName,'.noise'],'file') == 2))
        %%
        fprintf('Starting noise finding...\n');
        tic
        noise = RawDataNoiseEvaluationM(dataPath, saveFolder);
        
        fprintf('Time for noise evaluation %.2f seconds.\n', toc);
    else
        fprintf('Noise not requested or .noise file found - skipping raw data noise evaluation.\n');
    end
    
    
    %% Find spikes and make a .spikes file
    if tryToDo(2) &&...
            (force(2) || ~(exist([saveFolder,filesep,datasetName,'.spikes.mat'],'file') == 2))
        %%
        fprintf('Starting spike finding...\n');
        tic
        
        sigmaFileName = [saveFolder,filesep,datasetName,'.noise'];
        
        [spikes,ttlTimes,nSamples] = SpikeFindingM(dataPath, saveFolder, timeCommand, sigmaFileName);
        spikeSave = int32(spikes(:,1:2));
        save([saveFolder,filesep,datasetName,'.spikes.mat'],'spikeSave','ttlTimes','nSamples');
        %     save([saveFolder,filesep,datasetName,nameExt,'.spikes.mat'],'spikeSave','ttlTimes');
        
        fprintf('Time for spike finding %.2f seconds\n',toc);
    else
        fprintf('Spike not requested or .spikes.mat file found - skipping spike finding.\n');
    end
    
    
    %% Covariance calculation
    if tryToDo(3) &&...
            (force(3) || ~(exist([saveFolder,filesep,datasetName,'.cov.mat'],'file') == 2))
        %%
        fprintf('Starting covariance calculation...\n');
        tic
        
        if ~exist('spikeSave')
            load([saveFolder,filesep,datasetName,'.spikes.mat']);
        end
        
        [covMatrix,averages,totSpikes] = buildCovariances(double(spikeSave), dataPath, timeCommand);
        
        cfg = mVisionConfig(); covConfig = cfg.getCovConfig();
        if covConfig.whitening
            noiseSpikes = generateNoiseEvents(double(spikeSave));
           [noiseCovMatrix,~,~] = buildCovariances(noiseSpikes,dataPath,timeCommand);
           covMatrix = whitenMatrices(dataPath, totSpikes, covMatrix, noiseCovMatrix);
%          covMatrix = noiseCovMatrix;
        end
        save([saveFolder,filesep,datasetName,'.cov.mat'],'covMatrix','averages','totSpikes');
        
        fprintf('Time for covariance calculation %.2f seconds\n',toc);
    else
        fprintf('Cov not requested or .cov.mat file found - skipping covariance calculation.\n');
    end
    
    %% Eigenspikes Projections calculation
    if tryToDo(4) &&...
            (force(4) || ~(exist([saveFolder,filesep,datasetName,'.prj.mat'],'file') == 2))
        %%
        fprintf('Starting projections calculation...\n');
        tic
        if ~exist('covMatrix')
            load([saveFolder,filesep,datasetName,'.cov.mat']);
        end
        
        if ~exist('spikeSave')
            load([saveFolder,filesep,datasetName,'.spikes.mat']);
        end
        
        [projSpikes,eigenValues,eigenVectors,spikeTimes] = ...
            PCProj(dataPath, timeCommand, ...
            double(spikeSave), covMatrix, averages, totSpikes);
        
        save([saveFolder,filesep,datasetName,'.prj.mat'],'eigenValues','eigenVectors','spikeTimes','-v7.3');
        for el = 1:numel(projSpikes)
            eval(sprintf('projSpikes%u = projSpikes{el};',el));
            save([saveFolder,filesep,datasetName,'.prj.mat'],sprintf('projSpikes%u',el),'-append');
            eval(sprintf('projSpikes%u = [];',el));
            projSpikes{el} = []; % Progressive RAM clean-up
        end
        
        fprintf('Time for projections calculation %.2f seconds\n',toc);
    else
        fprintf('Prj not requested or .prj.mat file found - skipping projections calculation.\n');
    end
    
    
    %% Clustering
    if tryToDo(5) &&...
            (force(5) || ~(exist([saveFolder,filesep,datasetName,'.model.mat'],'file') == 2 &&...
            exist([saveFolder,filesep,datasetName,'.neurons.mat'],'file') == 2))
        %%
        fprintf('Starting clustering...\n')
        tic
        % No projections loader at this level
        % Loader of spikeTimes however
        if ~exist('spikeTimes','var')
            load([saveFolder,filesep,datasetName,'.prj.mat'],'spikeTimes');
        end
        % Clustering framework
        [clusterParams,neuronEls,neuronClusters,neuronSpikeTimes] =...
            PCClustering([saveFolder,filesep,datasetName,'.prj.mat'],spikeTimes);
        
        % Save a model file
        save([saveFolder,filesep,datasetName,'.model.mat'],'clusterParams');
        
        % Save a neurons-raw file (.neurons.mat)
        save([saveFolder,filesep,datasetName,'.neurons.mat'],'neuronEls','neuronClusters','neuronSpikeTimes');
        
        fprintf('Time for clustering %.2f seconds\n',toc);
    else
        fprintf('Clust not requested or .neurons|model.mat files found - skipping clustering.\n');
    end
    
    
    %% Cleaning and saving neurons in Vision compatible neuron file
    if tryToDo(6) &&...
            (force(6) || ~(exist([saveFolder,filesep,datasetName,'.neurons'],'file') == 2))
        %%
        fprintf('Removing duplicates and saving a vision-compatible .neurons file...\n')
        tic
        if ~exist('neuronSpikeTimes','var')
            load([saveFolder,filesep,datasetName,'.neurons.mat']);
        end
        
        [neuronEls, neuronClusters, neuronSpikeTimes] = ...
            duplicateRemoval(dataPath, saveFolder, datasetName, timeCommand, ...
            neuronEls, neuronClusters, neuronSpikeTimes);
        
        cfg = mVisionConfig(); dataConfig = cfg.getDataConfig();
        sampleRate = dataConfig.sampleRate;
        
        neuronSaver = NeuronSaverM(dataPath,saveFolder,datasetName,'',0);
        neuronSaver.pushAllNeurons(neuronEls, neuronClusters, neuronSpikeTimes);
        neuronSaver.close();
        
        fprintf('Neuron cleaning done.\n');
        
        fprintf('Time for cleaning and saving %.2f seconds\n', toc);
    else
        fprintf('Clean|Save not requested or .neurons file found - skipping cleaning|saving.\n');
    end
    
    %%
    fprintf('\n');
    fprintf('Total pipeline time %.2f seconds\n',toc(totalTime));
    fprintf([strrep(dataPath,'\','\\'),timeCommand,' finished\n']);
    fprintf('----------------------------------------------------------\n');
    
end
