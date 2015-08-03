%% Script file parallel caller
% Loads a list of datasets to process in ../fileList.input
% And processes them all in parallel on multiple cores
%%
config = mVisionConfig();
parConfig = config.getParConfig();

% parpool(parConfig.nWorkers);


%%
fileList = importdata(['..',filesep,'fileList.input']);
n = numel(fileList);

disp('parallelCaller: Entering file loop...');
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


%%
% delete(gcp);
exit;
