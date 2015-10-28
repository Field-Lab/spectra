function mVision(dataPath, saveFolder, timeCommand, movieXml, tryToDo, force)
    %MVISION Demonstrates workflow of mVision and forms the backbone
    % of serial neuron finding
    %
    % This script manages the successive steps of serial neuron finding
    % It can work in script mode or in function mode
    % Input:
    %   dataSetName: string in the format 'yyyy-mm-dd-x/dataxxx'
    %       The root folders for data and output are provided below
    %   timeCommand: vision style time command for selection within a dataset
    %       Such as (10-), (500-1000),(-100)
    %   mVision DOES NOT support dataset syntaxes such as data000(1000-)-data002(-100)
    %
    % Successively calls if requested || necessary the subroutines:
    % Noise evaluation - Spike Finding - Covariance Calculation
    % Projections Calculations - Clustering - Neuron Cleaning - Neuron File saving
    %
    % All intermediate files are stored in convenient .mat files, except:
    % .neurons-raw file saved under vision format and matlab format (.neurons.mat)
    % .neurons saved only in vision format
    %
    % VISION COMPATIBILITY
    % To run further analysis with vision, one should first call the calculations
    % "Generate globals file" then "Copy Raw Data Header To Globals"
    %
    % Author -- Vincent Deo -- Stanford University -- August 5, 2015
    
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
        % [filesep, 'concat'] is a special token here for concateated analysis
	throw(MException('','demoScript: data source foilder|file does not exist'));
    end
    
    mkdir(saveFolder);
    [~,datasetName,~] = fileparts(dataPath); % Catching dataset name as last part of saveFolder
    
    if exist(saveFolder,'file') == 2
        [~,datasetName,~] = fileparts(saveFolder);
    end
    
    % USER input - tryToDo -- won't do any task unless stated here
    nSteps = 7;
    % Steps in order:
    % --------- noise - spike - cov - prj - clust - save - stas ----------------------
    if isa(tryToDo,'char')
        tryToDo = lower(tryToDo);
        if strcmp(tryToDo,'all')
            tryToDo = ones(1,nSteps);
        end
        if strcmp(tryToDo,'nosta')
            tryToDo = [ones(1,nSteps-1),0];
        end
        if strcmp(tryToDo,'staonly')
            tryToDo = [ones(1,nSteps-1),0];
        end
        if strcmp(tryToDo,'noisetocov')
            tryToDo = [ones(1,3),zeros(1,nSteps-3)];
        end
        if strcmp(tryToDo,'prjtoneurons')
            tryToDo = [0,0,0,1,1,1,0];
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
        if isa(tryToDo,'char')
            throw(MException('','Invalid force argument'));
        end
    else
        validateattributes(force,{'numeric'},{'row','ncols',nSteps,'binary'},'','tryToDo',4);
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
        
        [neuronEls, neuronClusters, neuronSpikeTimes] =...
            duplicateRemoval(dataPath, saveFolder, datasetName, timeCommand, ...
            neuronEls, neuronClusters, neuronSpikeTimes);
        
        neuronSaver = NeuronSaverM(dataPath,saveFolder,datasetName,'');
        neuronSaver.pushAllNeurons(neuronEls, neuronClusters, neuronSpikeTimes);
        neuronSaver.close();
        
        fprintf('Neuron cleaning done.\n');
        
        fprintf('Time for cleaning and saving %.2f seconds', toc);
    else
        fprintf('Clean|Save not requested or .neurons file found - skipping cleaning|saving.\n');
    end
    
    %% Computing STAs
    if tryToDo(7) &&...
            (force(7) || ~(exist([saveFolder,filesep,datasetName,'.sta'],'file') == 2))
        %%
        fprintf('Computing STAs...\n')
        tic
        fprintf('STAs done.\n');
        f = filesep;
        if isunix
            c = ':'
        else
            c = ';'
        end
        system(['java -Xmx4g -Xss100m -cp ".',f,'vision',f,c,'.',f,'vision',f,'Vision.jar" edu.ucsc.neurobiology.vision.tasks.StaAnalysis ',...
            saveFolder,' ', movieXml, ' -c ..',f,'primate.xml']);
        
        fprintf('Time for STA computation %.2f seconds', toc);
    else
        fprintf('Clean|Save not requested or .neurons file found - skipping cleaning|saving.\n');
    end
    
    %%
    fprintf('\n');
    fprintf('Total pipeline time %.2f seconds\n',toc(totalTime));
    fprintf([strrep(dataPath,'\','\\'),timeCommand,' finished\n']);
    fprintf('----------------------------------------------------------\n');
    
end
