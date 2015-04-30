function covMatrix = buildCovariances(parameters, spikeFileName)
    % Build the covariance matrix for spikes around a given electrode
    % Input HashMap parameters should be the same given than for SpikeFindingM
    
    
    %% Imports
    import edu.ucsc.neurobiology.vision.electrodemap.*
    import edu.ucsc.neurobiology.vision.io.*
    import java.io.*
    
    %% Argument validation
    % Argument should be a java.util.HashMap<String,String> containing all relevant parameters for spike
    % finding
    validateattributes(parameters,{'java.util.HashMap'},{},'','parameters');
    p = parameters; % For concision
    
    %% Parsing and Storing Input HashMap
    rawDataSource = p.get('Raw_Data_Source'); % Actually at this point includes a command concatenated under the dataFileParser format: '.../data002(0-10)'
    sigmaPath = p.get('Sigma'); % .noise file
    outputPath = p.get('Analysis.Output_Path'); % Output path for the .spikes file
    
    meanTimeConstant = str2double(p.get('Mean Time Constant'));
    
    nLeftPoints = str2double(p.get('Analysis.Left Points'));
    nRightPoints = str2double(p.get('Analysis.Right Points'));
    minError = str2double(p.get('Analysis.Minimization Error'));
    spikeUse = p.get('Analysis.Spike To Use');
    
    %% New version: trying to get around using a MultipleCompressedSampleInputStream
    % Going to acquire samples through RawDataFile.getData()
    % Should allow to import larger number of samples and use array-based spikeFinding
    % TODO move that to a matlab file IO
    parser = DataFileStringParser(rawDataSource);
    datasets = parser.getDatasets();
    rawDataFile = RawDataFile(File(char(datasets(1))));
    startTimes = parser.getStartTimes();
    stopTimes = parser.getStopTimes();
    
    header = rawDataFile.getHeader();
    samplingRate = header.getSamplingFrequency();
    
    startSample = startTimes(1) * samplingRate;
    stopSample = stopTimes(1) * samplingRate;
    
    totalSamples = stopSample - startSample;
    
    
    
    
%% Java electrodemap setup if setElectrode is activated. Skipping this as set electrode is FALSE in vision config
    % if setElectrodes, array params from parameter HashMap
    % else array params from source header:
    % No idea so far what this is for
    packedArrayID = int32(header.getArrayID());
    binPackedArrayID = dec2bin(packedArrayID,32);
    arrayID = bin2dec(binPackedArrayID(17:32));
    arrayPart = bin2dec(binPackedArrayID(9:16));
    arrayNParts = bin2dec(binPackedArrayID(1:8));
    if arrayPart == 0
        arrayPart = 1;
        arrayNParts = 1;
    end
    
    electrodeMap = ElectrodeMapFactory.getElectrodeMap(packedArrayID);
    nElectrodes = electrodeMap.getNumberOfElectrodes();
    
    %% Creating spike file object
    spikeFile = SpikeFile(spikeFileName);
    spikeTimes = spikeFile.getSpikeTimes(); % returns all spike times in the file as a cell array
    % TODO should chunk as already 3MEG for 10 secs. Conclusion: .spikes are a waste
    
    
    
end