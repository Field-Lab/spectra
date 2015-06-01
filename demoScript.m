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
visionPath = [repoPath,'/vision/vision.jar'];
if ~exist('edu/ucsc/neurobiology/vision/io/RawDataFile','class')
    javaaddpath(visionPath)
end
javaaddpath('./vision');

% Get vision's config xml file
config = edu.ucsc.neurobiology.vision.Config([repoPath,'/vision/config.xml']);

% USER INPUT - Set up data and output folders
dataPath = 'X:\EJGroup_data\Data\2008-06-10-1\data000';
% dataPath = '/Volumes/Data/2013-04-30-3/data000'
timeCommand = '(0-10)';
saveFolder = 'X:\EJGroup_data\TestOut\2008-06-10-1\data000';
% saveFolder = '/home/vision/Vincent/mvision_outputs/2013-04-30-3/data000'

if ~(exist(dataPath,'file') == 2 || exist(dataPath,'file') == 7)
    throw(MException('demoScript','dataset folder does not exist'));
end
mkdir(saveFolder);
[~,dataSetName,~] = fileparts(dataPath); % Catching dataset name as last part of dataPath

%% Processing
%% Process noise and make a .noise file
if ~(exist([saveFolder,filesep,dataSetName,'.noise'],'file') == 2)
    disp('Starting noise finding...');
    tic
    noise = RawDataNoiseEvaluationM(dataPath,saveFolder);
    x = toc;
    disp(['Time for noise evaluation ', num2str(x), ' seconds.']);
else
    disp('.noise file found - skipping raw data noise evaluation.');
end

%% Find spikes and make a .spikes file
if ~(exist([saveFolder,filesep,dataSetName,'.spikes'],'file') == 2)
    disp('Starting spike finding...');
    tic
    
    sigmaFileName = [saveFolder,filesep,dataSetName,'.noise'];
    parameters = spikeFindingSetup([dataPath,timeCommand],saveFolder,sigmaFileName,config);
    
    SpikeFindingM(parameters);
    
    x = toc;
    disp(['Time for spike finding ', num2str(x), ' seconds']);
else
    disp('.spikes file found - skipping spike finding.');
end

%% Covariance calculation
if ~(exist([saveFolder,filesep,dataSetName,'.cov.mat'],'file') == 2)
    disp('Starting covariance calculation...');
    tic
    
    sigmaFileName = [saveFolder,filesep,dataSetName,'.noise'];
    parameters = spikeFindingSetup([dataPath,timeCommand],saveFolder,sigmaFileName,config);
    spikeFileName = [saveFolder,filesep,dataSetName,'.spikes'];
    
    [covMatrix,averages,totSpikes] = buildCovariances(parameters, spikeFileName);
    
    save([saveFolder,filesep,dataSetName,'.cov.mat'],'covMatrix','averages','totSpikes');
    x = toc;
    disp(['Time for covariance calculation ', num2str(x), ' seconds']);
else
    disp('.cov.mat file found - skipping covariance calculation.');
end

%% Eigenspikes Projections calculation
if ~(exist([saveFolder,filesep,dataSetName,'.prj.mat'],'file') == 2)
    disp('Starting projections calculation...');
    tic
    if ~exist('covMatrix')
        load([saveFolder,filesep,dataSetName,'.cov.mat']);
    end
    
    sigmaFileName = [saveFolder,filesep,dataSetName,'.noise'];
    spikeFileName = [saveFolder,filesep,dataSetName,'.spikes'];
    parameters = spikeFindingSetup([dataPath,timeCommand],saveFolder,sigmaFileName,config);
    
    [projSpikes,eigenValues,eigenVectors] = ...
        PCProj(parameters, spikeFileName, covMatrix, averages, totSpikes, 10);
    
    save([saveFolder,filesep,dataSetName,'.prj.mat'],'projSpikes','eigenValues','eigenVectors');
    
    x = toc;
    disp(['Time for projections calculation ', num2str(x), ' seconds']);
else
    disp('.prj.mat file found - skipping projections calculation.');
end