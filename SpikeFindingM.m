function [ ] = SpikeFindingM( parameters )
%SPIKEFINDINGM Matlab implementation of the Spike Finding algorithm
%   Takes on after noise finding has been done (or not)
% Computes the spikesproperties over the electrode array and stores
% them in a .spike file
% In a format matching the current behavior of vision.

import edu.ucsc.neurobiology.vision.electrodemap.*

%% Argument validation
% Argument should be a java.util.HashMap<String,String> containing all relevant parameters for spike
% finding
validateattributes(parameters,{'java.util.HashMap'},{},'','parameters');
p = parameters; % For concision

%% Parsing and Storing input HashMap
rawDataSource = p.get('Raw_Data_Source'); % Actually at this point includes a command concatenated under the dataFileParser format: '.../data002(0-10)'
sigmaPath = p.get('Sigma'); % .noise file

spikeThreshold = str2double(p.get('Spike Threshold'));
ttlThreshold = str2double(p.get('TTL Threshold'));
meanTimeConstant = str2double(p.get('Mean Time Constant'));

bufferSizeInBytes = 1024*str2double(p.get('Buffer Size (Kb)')); % NOT casting in int, because Matlab deals poorly with int divisions...
bufferSizeInBytes = round(bufferSizeInBytes/770)*770; % One sample in buffer takes 770 bytes - reducing buffer to int number of samples. 770 is 513 max electrodes * 12 bits number

nBuffers = str2double(p.get('Buffers Count')); % id line above

% Note : no option waitForData accepted for now. Input can't be a live feed or a being-written file.

% Instantiate data feed
% sampleInputStream = MultipleCompressedSampleInputStreamM(rawDataSource, bufferSizeInBytes, nBuffers);

% Extracting header
% header = sampleInputStream.getJavaHeader();
% nSamples = header.getNumberOfSamples();

%% New version: trying to get around using a MultipleCompressedSampleInputStream
% Going to acquire samples through RawDataFile.getData()
% Should allow to import larger number of samples and use array-based spikeFinding
% TODO move that to a matlab file IO
parser = DataFileStringParser(dataPath);
datasets = parser.getDatasets();
rawDataFile = RawDataFile(File(char(datasets(1))));
startTimes = parser.getStartTimes();
startTime = startTimes(1);
stopTimes = parser.getStopTimes();
stopTime = stopTimes(1);

header = rawDataFile.getHeader();
samplingRate = 20000;
totalSamples = (stopTime - startTime) * samplingRate;

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

%% Get sigmas
% sigma = getSigmas(sigmaPath, nElectrodes)
% Implemeting functionality of SpikeFinding.getSigmas(String fileNameOrValue, nElectrodes);
% TODO Check better and throw an error
if isempty(str2num(sigmaPath)) % sigmaPath is a path string
    sigma = load(sigmaPath,'-ascii');
    validateattributes(sigma,{'numeric'},{'ncols',1,'nrows',nElectrodes},'','loaded sigma file');
    % File path error
else % sigmaPath is a value
    sigma = sigmaPath * ones(nElectrodes,1);
end
sigma(1) = 100; % TTL electrode
sigma = sigma * spikeThreshold;


%% Raw Data Saving
% Skipped - this version meant to work on already recorded data

%% Create the Spiker Finder and heirs
spikeFinderM = SpikeFinderM(electrodeMap, sigma, ttlThreshold, meanTimeConstant);

%% Analysis
% Skipped - this version meant to do spike saving

%% Diagnostic - None
%% GUI Stuff - None
%% Get the machine ready
%% Spike Finding
% Covers the SpikeFinder and SpikeBuffer in Matlab.
% Then passes spikes down the SpikeSaver


end

