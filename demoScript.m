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

% USER INPUT - Path to vision jar - add to java path
visionPath = 'C:/Users/Vincent/Documents/EJGroup/Java/vision7/Vision.jar';
if ~exist('edu/ucsc/neurobiology/vision/io/RawDataFile','class')
    javaaddpath(visionPath)
end

% Get vision's config xml file
[visionFolder,~,~] = fileparts(visionPath);
config = edu.ucsc.neurobiology.vision.Config([visionFolder,filesep,'config.xml']);

% USER INPUT - Set up data and output folders
dataPath = 'X:\EJGroup_data\Data\2005-04-26-0\data002';
saveFolder = 'X:\EJGroup_data\TestOut\2005-04-26-0\data002';
if ~(exist(dataPath,'file') == 2 || exist(dataPath,'file') == 7)
    throw(MException('demoScript','dataset folder does not exist'));
end
mkdir(saveFolder);
[~,dataSetName,~] = fileparts(dataPath); % Catching dataset name as last part of dataPath

%% Processing
%% Process noise and make a .noise file
if ~(exist([saveFolder,filesep,dataSetName,'.noise'],'file') == 2)
    noise = RawDataNoiseEvaluationM(dataPath,saveFolder);
end

%% Find spikes and make a .spikes file
if ~(exist([saveFolder,filesep,dataSetName,'.spikes'],'file') == 2)
    sigmaFileName = [saveFolder,filesep,dataSetName,'.noise'];
    parameters = spikeFindingSetup([dataPath,'(0-10)'],saveFolder,sigmaFileName,config);
    SpikeFindingM(parameters);
end

