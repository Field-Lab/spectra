function outputHashMap = spikeFindingSetup(rawDataFile,outputPath,sigmaFileName,visionConfig)
    % SPIKEFINDINGSETUP Utilitarian to generate parameter HashMap
    % for vision SpikeFinding given a vision Config object
    
    import edu.ucsc.neurobiology.vision.*
    import edu.ucsc.neurobiology.vision.anf.*
    import java.util.*
    import java.lang.String
    
    % Sigma argument ?
    % Refer to edu.ucsc.neurobiology.vision.tasks.RunScript.
    % spikeFinding(String, String, String, AnalysisToDo,boolean,Config)
    
    % Parameter check
    validateattributes(visionConfig,{'edu.ucsc.neurobiology.vision.Config'},{},'','visionConfig');
    % TODO
    if false % Path validity for rawDataFile
        ME = MException('SpikeFindingSetup:IllegalArgumentException',...
            ['argument rawDataFile ', rawDataFile,' is not a valid path.']);
        throw(ME);
    end
    % TODO
    if false % Path validity for outputPath -- Not implemented -- See exist/ check dataFileNameParser for specs.
        ME = MException('SpikeFindingSetup:IllegalArgumentException',...
            ['argument outputPath ',outputPath,' is not a valid path.']);
        throw(ME);
    end
    % TODO
    if false % Path validity for sigmaFileName -- Same as above
        ME = MException('SpikeFindingSetup:IllegalArgumentException',...
            ['argument sigmaFilename ',sigmaFileName,' is not a valid path.']);
        throw(ME);
    end
    
    %% Compatibilized java code
    
    cName = 'Spike Finding';
    d = visionConfig.getParameterList(cName);
    
    outputHashMap = HashMap();
    outputHashMap.put('Raw_Data_Source', rawDataFile);
    outputHashMap.put('Buffer Size (Kb)', d.get('Buffer Size (Kb)'));
    outputHashMap.put('Buffers Count', d.get('Buffers Count'));
    outputHashMap.put('Spike Threshold', d.get('Spike Threshold'));
    outputHashMap.put('Sigma', sigmaFileName);
    outputHashMap.put('TTL Threshold', d.get('TTL Threshold'));
    outputHashMap.put('Mean Time Constant', d.get('Mean Time Constant'));
    if d.containsKey('waitForData')
        outputHashMap.put('waitForData', d.get('waitForData'));
    end
    
    outputHashMap.put('Diagnostic Plots', 'false');
    outputHashMap.put('saveRawData', 'false');
    
    outputHashMap.put('Set Electrodes', d.get('Set Electrodes'));
    outputHashMap.put('Set Electrodes.arrayID', d.get('Set Electrodes.arrayID'));
    outputHashMap.put('Set Electrodes.arrayPart', d.get('Set Electrodes.arrayPart'));
    outputHashMap.put('Set Electrodes.arrayNParts', d.get('Set Electrodes.arrayNParts'));
    outputHashMap.put('Set Electrodes.flipX', d.get('Set Electrodes.flipX'));
    outputHashMap.put('Set Electrodes.flipY', d.get('Set Electrodes.flipY'));
    
    % Case analysisToDo = SAVE_SPIKES in java
    % No interest in doing nothing
    % or doing covariance analysis concurrent to spike finding for now.
    % Below is a trick to recover a nested enum value:
    analysisToDo = javaMethod('valueOf', 'edu.ucsc.neurobiology.vision.anf.SpikeFinding$AnalysisToDo', 'SAVE_SPIKES');
    
    outputHashMap.put('Analysis', 'true');
    outputHashMap.put('Analysis.Analysis To Do', String(num2str(analysisToDo.ordinal())));
    outputHashMap.put('Analysis.Left Points', String(d.get('Analysis.Left Points')));
    outputHashMap.put('Analysis.Right Points', String(d.get('Analysis.Right Points')));
    outputHashMap.put('Analysis.Minimization Error', d.get('Analysis.Minimization Error'));
    outputHashMap.put('Analysis.Spike To Use', d.get('Analysis.Spike To Use'));
    % if (!unWhitenedCovariances) {
    %    outputHashMap.put('Analysis.Minimum Noise Events', d.get('Analysis.Minimum Noise Events'));
    %    } else {
    outputHashMap.put('Analysis.Minimum Noise Events', String('0'));
    %        }    
    outputHashMap.put('Analysis.Electrode Usage', String(num2str(d.get('Analysis.Electrode Usage'))));
    outputHashMap.put('Analysis.Output_Path', outputPath);
    
    
end