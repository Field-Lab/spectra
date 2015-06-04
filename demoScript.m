%%%%%%% Demo Script %%%%%%%%%
% This script demonstrates the data flow in mVision
% Starting at raw data files and ending at .neurons file.
%
% Work in progress
% Current state: clustering not implemented
% writing dummy clusters into .neurons
%
% Vincent Deo - Stanford University - 05/29/2015


%% SETUP
clear;
% Add subfolders to the matlab path
addpath(genpath('./'));

% generate repository path
repoPath = pwd;

% USER INPUT - Path to vision jar - add to java path
visionPath = [repoPath,'/vision/Vision.jar'];
if ~exist('edu/ucsc/neurobiology/vision/io/RawDataFile','class')
    javaaddpath(visionPath)
end
javaaddpath('./vision');

% Get vision's config xml file
config = edu.ucsc.neurobiology.vision.Config([repoPath,'/vision/config.xml']);

% USER INPUT - Set up data and output folders
% dataPath = 'X:\EJGroup_data\Data\2008-06-10-1\data000'
dataPath = '/Volumes/Data/2013-04-30-3/data000'
timeCommand = '(0-10)';
% saveFolder = 'X:\EJGroup_data\TestOut\2008-06-10-1\data000Matlab'
saveFolder = '/home/vision/Vincent/mvision_outputs/2013-04-30-3/data000'

% USER input - FORCE rewriting output even if files are found
force = 0;
% 0 force all - 1 force from spikes - 2 force from cov - 3 force from proj
% 4 force from clustering and cleaning - 5 force vision .neuron rewrite
% 6 force none



if ~(exist(dataPath,'file') == 2 || exist(dataPath,'file') == 7)
    throw(MException('demoScript','dataset folder does not exist'));
end
mkdir(saveFolder);
[~,datasetName,~] = fileparts(dataPath); % Catching dataset name as last part of dataPath


%% Process noise and make a .noise file
if force <= 0 || ~(exist([saveFolder,filesep,datasetName,'.noise'],'file') == 2)
    %%
    disp('Starting noise finding...');
    tic
    noise = RawDataNoiseEvaluationM(dataPath,saveFolder);
    x = toc;
    disp(['Time for noise evaluation ', num2str(x), ' seconds.']);
else
    disp('.noise file found - skipping raw data noise evaluation.');
end


%% Find spikes and make a .spikes file
if force <= 1 || ~(exist([saveFolder,filesep,datasetName,'.spikes'],'file') == 2)
    %%
    disp('Starting spike finding...');
    tic
    
    sigmaFileName = [saveFolder,filesep,datasetName,'.noise'];
    parameters = spikeFindingSetup([dataPath,timeCommand],saveFolder,sigmaFileName,config);
    
    SpikeFindingM(parameters);
    
    x = toc;
    disp(['Time for spike finding ', num2str(x), ' seconds']);
else
    disp('.spikes file found - skipping spike finding.');
end


%% Covariance calculation
if force <= 2 || ~(exist([saveFolder,filesep,datasetName,'.cov.mat'],'file') == 2)
    %%
    disp('Starting covariance calculation...');
    tic
    
    sigmaFileName = [saveFolder,filesep,datasetName,'.noise'];
    parameters = spikeFindingSetup([dataPath,timeCommand],saveFolder,sigmaFileName,config);
    spikeFileName = [saveFolder,filesep,datasetName,'.spikes'];
    
    [covMatrix,averages,totSpikes] = buildCovariances(parameters, spikeFileName);
    
    save([saveFolder,filesep,datasetName,'.cov.mat'],'covMatrix','averages','totSpikes');
    x = toc;
    disp(['Time for covariance calculation ', num2str(x), ' seconds']);
else
    disp('.cov.mat file found - skipping covariance calculation.');
end


%% Eigenspikes Projections calculation
if force <= 3 || ~(exist([saveFolder,filesep,datasetName,'.prj.mat'],'file') == 2)
    %%
    disp('Starting projections calculation...');
    tic
    if ~exist('covMatrix')
        load([saveFolder,filesep,datasetName,'.cov.mat']);
    end
    
    sigmaFileName = [saveFolder,filesep,datasetName,'.noise'];
    spikeFileName = [saveFolder,filesep,datasetName,'.spikes'];
    parameters = spikeFindingSetup([dataPath,timeCommand],saveFolder,sigmaFileName,config);
    
    [projSpikes,eigenValues,eigenVectors,spikeTimes] = ...
        PCProj(parameters, spikeFileName, covMatrix, averages, totSpikes, 10);
    
    save([saveFolder,filesep,datasetName,'.prj.mat'],'projSpikes','eigenValues','eigenVectors','spikeTimes');
    
    x = toc;
    disp(['Time for projections calculation ', num2str(x), ' seconds']);
else
    disp('.prj.mat file found - skipping projections calculation.');
end


%% Clustering + Neuron Cleaning
if force <= 4 || ~(exist([saveFolder,filesep,datasetName,'.model.mat'],'file') == 2 &&...
        exist([saveFolder,filesep,datasetName,'.neurons.mat'],'file') == 2)
    %%
    disp('Starting clustering and neuron cleaning...')
    tic
    if ~exist('projSpikes','var')
        load([saveFolder,filesep,datasetName,'.prj.mat']);
    end
    % Do clustering stuff
    
    % Write a .model.mat containing the clustering information
    
    % Separate the neurons in format [neuronID, neuronSpikeTimes]
    % Do some neuron cleaning if wanted
 
    [clusterParams,neuronEls,neuronClusters,neuronSpikeTimes] = PCClustering(projSpikes, spikeTimes);
    
    save([saveFolder,filesep,datasetName,'.model.mat'],'clusterParams');
    
%     neuronEls = [1;1;3;4];
%     neuronClusters = [1;2;1;1];
%     neuronSpikeTimes = {[1,2,3,4];[11,12,13,14];[21,22,23,24];[25,26,27,29]};
    save([saveFolder,filesep,datasetName,'.neurons.mat'],'neuronEls','neuronClusters','neuronSpikeTimes');
    
    x = toc;
    disp(['Time for clustering ', num2str(x), ' seconds']);
else
    disp('.prj.mat file found - skipping projections calculation.');
end


%% Saving neurons in Vision compatible neuron file
if force <=5 || ~(exist([saveFolder,filesep,datasetName,'.neurons'],'file') == 2)
    %%
    disp('Saving a vision-compatible .neurons file...')
    tic
    if ~exist('neuronSpikeTimes','var')
        load([saveFolder,filesep,datasetName,'.neurons.mat']);
    end
    neuronSaver = NeuronSaverM(dataPath,saveFolder,datasetName);
    for i = 1:numel(neuronEls)
        el = neuronEls(i);
        neuronSaver.addNeuron(el,...
            neuronSaver.getNeuronID(el,neuronClusters(i)),...
            neuronSpikeTimes{i});
    end
    
    x = toc;
    disp(['Time for saving ', num2str(x), ' seconds']);
else
    disp('.neurons file found - skipping saving.');
end


