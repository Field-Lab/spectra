function CC = corrValFunction( dataFolder, datasetName )
    %CORRVALFUNCTION Summary of this function goes here
    %   Detailed explanation goes here
    neurPath = [dataFolder,filesep,datasetName,'.neurons.mat'];
    corrPath = [dataFolder,filesep,datasetName,'.corr.mat'];
    load(neurPath);
    CC = eye(numel(neuronSpikeTimes));
    for i = 1:numel(neuronSpikeTimes)
        for j = (i+1):numel(neuronSpikeTimes)
            CC(i,j) = edu.ucsc.neurobiology.vision.anf.NeuronCleaning.getCorrVal(neuronSpikeTimes{i},neuronSpikeTimes{j},10);
            CC(j,i) = CC(i,j);
        end
    end
    save(corrPath,'CC');
end

