% Script comparatorScript
% Loads a list of datasets to process in ../fileList.input
% Sets up folders and compares reference raw neurons vs new neurons
%
% 
% Author -- Vincent Deo -- Stanford University -- August 27, 2015

%% Load config - Start parallel pool
config = mVisionConfig();
parConfig = config.getParConfig();

% parpool(parConfig.nWorkers);

%% Load input
fileList = textread(['..',filesep,'fileList.input'],'%s');
n = numel(fileList);
scores = cell(n,1);

disp('Neuron comparator: Entering file loop...');
% Start parallel loop on files
% parfor
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
%% Store output in table format and save
scoreCat = horzcat(scores{:});

m = [cell2mat(cellfun(@(X) mean(X),scores,'uniformoutput',false));mean(scoreCat)];
s = [cell2mat(cellfun(@(X) std(X),scores,'uniformoutput',false));std(scoreCat)];
m2 = [cell2mat(cellfun(@(X) mean(X(X ~= 0)),scores,'uniformoutput',false));mean(scoreCat(scoreCat ~= 0))];
s2 = [cell2mat(cellfun(@(X) std(X(X ~= 0)),scores,'uniformoutput',false));std(scoreCat(scoreCat ~= 0))];

resultTable = ...
    table([fileList;'Global'],m,s,m2,s2,...
    'VariableNames',{'Dataset','avScores','std','avScoresNoZeros','stdNoZeros'})

save('/Volumes/Lab/Projects/spikesorting/mvision/scores/scoresGaussianClustering.mat','scores','resultTable','-v7.3');

%% Delete parallel pool, quit.
% delete(gcp);
exit;
