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
config = mVisionConfig();
parConfig = config.getParConfig();

parpool(parConfig.nWorkers);

%%
fileList = textread(['..',filesep,'fileList.input'],'%s');
n = numel(fileList);

% parfor
for k = 1:n
    try
        demoScript(fileList{k},'');
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
exit;
