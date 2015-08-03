%% Script comparatorScript
% Loads a list of datasets to process in ../fileList.input
% Sets up folders and compares reference raw neurons vs new neurons
%%
config = mVisionConfig();
parConfig = config.getParConfig();

% parpool(parConfig.nWorkers);

%%
fileList = importdata(['..',filesep,'fileList.input']);
n = numel(fileList);
scores = cell(n,1);

disp('Neuron comparator: Entering file loop...');
% for
for k = 1:n
    try
        scores{k} = neuronComparator(fileList{k});
    catch error
        disp(['Error in file ',fileList{k}]);
        disp(error);
        for i = 1:numel(error.stack)
            disp(error.stack(i));
        end
    end
end

save('/Volumes/Lab/Projects/spikesorting/mvision/scores/scoresSpectralClustering.mat',scores,'-v7.3');

%%
% delete(gcp);
exit;
