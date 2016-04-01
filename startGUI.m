function startGUI(path)
    % STARTGUI starter for cluster viewer GUI
    % Input:
    %       path: path to analysis folder
    %
    % Handles java path
    % Passes the path argument down to clusterEditGUI
    % In case of non-existing classes errors, or "java objects exist" warnings
    % run a "clear java" and try again.
    
    %% MATLAB VERSION CHECK
    mver = ver('MATLAB');
    if str2double(mver.Version) < 8.5
        throw(MException('','startGUI.m: Spectra requires MATLAB R2015a or later.'));
    end
    
    %% Start the GUI
    addpath(genpath(['.',filesep]));
    if ~exist('edu.ucsc.neurobiology.vision.Vision','class')
        javaaddpath ./vision/Vision.jar -end
        javaaddpath ./vision/
        javaaddpath ./duplicateRemoval/java_EI_comparison/
        javaaddpath ./clusterUtil/
    end
    
    ClusterEditGUI(path);
end
