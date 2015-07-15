% parpool(4);

%%
fileList = importdata(['..',filesep,'fileList.input']);
n = numel(fileList);

parfor k = 1:n
    neuronComparator(fileList{k});
end


%%
% delete(gcp);
% exit;
