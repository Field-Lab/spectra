%%
parpool(4);

%%
fileList = importdata(['..',filesep,'fileList.input']);
n = numel(fileList);

parfor k = 1:n
    disp(['Starting on ',fileList{k}]);
    pause(0.1);
    neuronComparator(fileList{k});
end


%%
delete(gcp);
exit;
