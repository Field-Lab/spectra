%% Script file parallel caller
%
% Loads a list of datasets to process in ../fileList.input
% And processes them all in sequential or parallel
%
% Warning: A parfor should not be used here if another (more efficient) is used
% in subroutines (e.g. in PCClustering).
%
% Author -- Vincent Deo -- Stanford University -- August 5, 2015

%% Load mVision configuration - start parallel pool if needed
CONFIG_TAG = ''; % Put your config tag here. It is passed down to mVision
% But not globalized from this level.
addpath('./cfg/');
cfg = mVisionConfig(CONFIG_TAG);
parConfig = cfg.getParConfig();

parpool(parConfig.nWorkers);

%%
fileList = textread(['/Volumes/Lab/Projects/spikesorting/mvision',filesep,'fileList2.input'],'%s');
movieList = textread(['/Volumes/Lab/Projects/spikesorting/mvision',filesep,'movies2.input'],'%s');
n = numel(fileList);

prefix = '/Volumes/Archive/'
outputPrefix = '/Volumes/Lab/Projects/spikesorting/mvision/outputsFeb16/'
moviePrefix = '/Volumes/Analysis/stimuli/white-noise-xml/'

% parfor
parfor k = 1:n
% for k = 1:1
%    try
        mVision([prefix,fileList{k}],... % Data folder
            [outputPrefix,fileList{k}],... % Output folder
            '',... % Time tag
            [moviePrefix, movieList{k}, '.xml'],... % movieXML
            [0 0 0 0 0 1 1 1],... % Requested computation
            'none',... % force overwrite
            CONFIG_TAG); % Config tag
        
%    catch error
%       disp(['Error in file ',fileList{k}]);
%       disp(error);
%       for i = 1:numel(error.stack)
%           disp(error.stack(i));
%       end
%    end
end

%% Clean and close
delete(gcp);

%% profiler
% profile off;
% p = profile('info');
% save([outputPrefix,'profile'],'p');

exit;
