function CC = corrValFunction( dataFolder, datasetName )
    %CORRVALFUNCTION Summary of this function goes here
    %   Detailed explanation goes here
    neurPath = [dataFolder,filesep,datasetName,'.neurons'];
    corrPath = [dataFolder,filesep,datasetName,'.corr.mat'];
    nrf = edu.ucsc.neurobiology.vision.io.NeuronFile(neurPath);
    IDList = nrf.getIDList();
    
    CC = eye(numel(IDList));
    for i = 1:numel(IDList)
        for j = (i+1):numel(IDList)
            CC(i,j) = edu.ucsc.neurobiology.vision.anf.NeuronCleaning.getCorrVal(nrf,IDList(i),IDList(j),10);
            CC(j,i) = CC(i,j);
        end
    end
    save(corrPath,'CC');
    nrf.close();
end

