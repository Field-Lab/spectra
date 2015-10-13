function demoScript(varargin)
    %DEMOSCRIPT Demonstrates workflow of mVision and forms the backbone
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
    
    if ~(nargin == 0 || nargin == 2)
        throw(MException('','Need 0 variables if script or 2 variables (''yyyy-mm-dd-x/data00x'',timeCommand) if function mode'));
    end
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
    if nargin ~= 2 % Script mode - put your dataset to process
        dataPath = 'X:\EJGroup_data\Data\2005-04-26-0\data002'
        % dataPath = '/Volumes/Data/2013-04-30-3/data001'
        
        timeCommand = ''
        
        saveFolder = 'X:\EJGroup_data\sshed_Results_stas_eis\copyForTest-newCleaning\2005-04-26-0\data002'
        % saveFolder = '/home/vision/vincent/outputs/2013-04-30-3/data001'
    else % function mode
        % USER - Sethere root folder for data/output of your list of datasets
        dataPath = ['/Volumes/Archive/',varargin{1}]
        saveFolder = ['/Volumes/Lab/Projects/spikesorting/mvision/outputsSpectral/',varargin{1}]
        timeCommand = varargin{2}
    end
    
    % DEBUG - additional saved file dataset name extension
    nameExt = '';
    
    % USER input - tryToDo -- won't do any task unless stated here
    % --------- noise - spike - cov - prj - clust - save ----------------------
    tryToDo =  [  0   ,   0   ,  0  ,  0  ,   0   ,   1  ];
    % USER input - force -- rewriting output even if files are found
    % --------- noise - spike - cov - prj - clust - save ----------------------
    force =    [  0   ,   0   ,  0  ,  0  ,   0   ,   1  ];
    
    if ~(exist(dataPath,'file') == 2 || exist(dataPath,'file') == 7)
        throw(MException('','demoScript: data source folder|file does not exist'));
    end
    mkdir(saveFolder);
    [~,datasetName,~] = fileparts(dataPath); % Catching dataset name as last part of dataPath
    if exist(dataPath,'file') == 2
        datasetName = datasetName(1:7);
    end
    
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
        
        [spikes,ttlTimes] = SpikeFindingM(dataPath, saveFolder, timeCommand, sigmaFileName);
        spikeSave = int32(spikes(:,1:2));
        save([saveFolder,filesep,datasetName,'.spikes.mat'],'spikeSave','ttlTimes');
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
    
    %%
    fprintf('\n');
    fprintf('Total pipeline time %.2f seconds\n',toc(totalTime));
    fprintf([strrep(dataPath,'\','\\'),timeCommand,' finished\n']);
    fprintf('----------------------------------------------------------\n');
    
end
