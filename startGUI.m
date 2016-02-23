function startGUI(path)
    
    addpath(genpath(['.',filesep]));
    % Add warning filters
    javaaddpath ./vision/Vision.jar -end
    % For custom use only, if using java changes in the local path
    % Uncompiled yet in the jar
    javaaddpath ./vision/
    javaaddpath ./duplicateRemoval/java_EI_comparison/
    javaaddpath ./clusterUtil/
    
    ClusterEditGUI(path);
end