cd /Volumes/Lab/Projects/spikesorting/mvision/mvision
addpath(genpath(['.',filesep]));

% Add warning filters
javaaddpath ./vision/Vision.jar -end
javaaddpath ./vision/
javaaddpath ./duplicateRemoval/java_EI_comparison/
javaaddpath ./clusterUtil/
javaclasspath('-dynamic')

ClusterEditGUI('/Volumes/Lab/Projects/spikesorting/mvision/whiteCompTest/withNewCleaning/data002');