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
profile on;

config = mVisionConfig();
parConfig = config.getParConfig();

parpool(parConfig.nWorkers);

%%
fileList = textread(['..',filesep,'fileList.input'],'%s');
movieList = textread(['..',filesep,'movies.input'],'%s');
n = numel(fileList);

prefix = '/Volumes/Archive/'
outputPrefix = '/Volumes/Lab/Projects/spikesorting/mvision/outputsSpectral/'
movieprefix = '/Volumes/Analysis/stimuli/white-noise-xml/'

% parfor
% for k = 1:n
for k = 1:1
    try
        mVision([prefix,fileList{k}],[outputPrefix,fileList{k}],'',[movieprefix,movieList{k}],'all','all');
    catch error
       disp(['Error in file ',fileList{k}]);
       disp(error);
       for i = 1:numel(error.stack)
           disp(error.stack(i));
       end
    end
end

%% Clean and close
delete(gcp);

%% profiler
profile off;
p = profile('info');
save([outputPrefix,'profile'],'p');

exit;
