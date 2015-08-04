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
%%
scoreCat = horzcat(scores{:});

m = [cell2mat(cellfun(@(X) mean(X),scores,'uniformoutput',false));mean(scoreCat)];
s = [cell2mat(cellfun(@(X) std(X),scores,'uniformoutput',false));std(scoreCat)];
m2 = [cell2mat(cellfun(@(X) mean(X(X ~= 0)),scores,'uniformoutput',false));mean(scoreCat(scoreCat ~= 0))];
s2 = [cell2mat(cellfun(@(X) std(X(X ~= 0)),scores,'uniformoutput',false));std(scoreCat(scoreCat ~= 0))];

resultTable = ...
    table([fileList;'Global'],m,s,m2,s2,...
    'VariableNames',{'Dataset','avScores','std','avScoresNoZeros','stdNoZeros'})
%%
save('/Volumes/Lab/Projects/spikesorting/mvision/scores/scoresGaussianClustering.mat','scores','resultTable','-v7.3');

%%
% delete(gcp);
exit;
