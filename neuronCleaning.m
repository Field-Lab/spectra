function neuronCleaning( neuronFileName )
    %NEURONCLEANING Calls vision's neuron cleaning
    %   Loads an existing neuron file
    %   Cleans neurons and overwrites.
    
    p = java.util.HashMap();
    p.put('Neuron_File',neuronFileName);
    
    config = mVisionConfig();
    cleanConfig = config.getCleanConfig();
    
    p.put('Minimun Number of Spikes',num2str(cleanConfig.minSpikes));
    p.put('Maximum Contamination',num2str(cleanConfig.maxCont));
    p.put('Coincidence Time',num2str(cleanConfig.coincTime));
    p.put('Maximum Correlation',num2str(cleanConfig.maxCorr));
    
    cleaningCalculation = edu.ucsc.neurobiology.vision.anf.NeuronCleaning();
    cleaningCalculation.setParameters(p);
    cleaningCalculation.startCalculation();
end

