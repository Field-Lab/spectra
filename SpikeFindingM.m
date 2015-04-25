function [ ] = SpikeFindingM( parameters )
%SPIKEFINDINGM Matlab implementation of the Spike Finding algorithm
%   Takes on after noise finding has been done (or not)
% Computes the spikesproperties over the electrode array and stores
% them in a .spike file
% In a format matching the current behavior of vision. %% Note: offset of 20000 samples, probably
% due to initialization of spikeFinder.

import edu.ucsc.neurobiology.vision.electrodemap.*
import edu.ucsc.neurobiology.vision.io.*
import java.io.*

%% Argument validation
% Argument should be a java.util.HashMap<String,String> containing all relevant parameters for spike
% finding
validateattributes(parameters,{'java.util.HashMap'},{},'','parameters');
p = parameters; % For concision

%% Parsing and Storing input HashMap
rawDataSource = p.get('Raw_Data_Source'); % Actually at this point includes a command concatenated under the dataFileParser format: '.../data002(0-10)'
sigmaPath = p.get('Sigma'); % .noise file
outputPath = p.get('Analysis.Output_Path'); % Output path for the .spikes file

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
spikeBufferM = SpikeBufferM();
spikeSaverM  = SpikeSaverM(header, outputPath);

%% Analysis
% Skipped - this version meant to do spike saving

%% Diagnostic - None
%% GUI Stuff - None
%% Get the machine ready
%% Spike Finding

lastSampleLoaded = 0; % Index of last sample read, initialize @ 0
currentSample = 0;

if totalSamples >= samplingRate % We need 1 sec buffer to initialize Spikefinder with
    rawData = rawDataFile.getData(startSample, samplingRate)'; % Get a read of 1 second
    lastSampleLoaded = samplingRate; % Do not update? - Seems that Vision starts at sample 1 even after initializing
    
    spikeFinderM.initialize(rawData);
else
   throw(MException('SpikeFindingM',...
                    'Total number of samples is insufficient to initialize SpikeFinder')); 
end

while lastSampleLoaded < stopSample-1 % stopSample should be the first sample not loaded
   rawData = rawDataFile.getData(lastSampleLoaded+1, min(samplingRate,stopSample-lastSampleLoaded-1))';
   lastSampleLoaded = lastSampleLoaded + min(samplingRate,stopSample-lastSampleLoaded-1);
   
   for i = 1:size(rawData,2)
       currentSample = currentSample + 1;
       spikeBufferM.addSpikes(spikeFinderM.processSample(double(rawData(:,i))));
       s = spikeBufferM.getSpikes(currentSample);
       spikeSaverM.processSpikes(s);
   end
end

spikeSaverM.processSpikes(spikeBufferM.getAllSpikes());
spikeSaverM.finishSpikeProcessing();
% Initialize stuff with this read


% Loop loading buffers
% Pass buffers to spikeFinder
% Pass Spikes to spikes buffer
% Pass Spikes to spikeSaver

% Get All remaining Spikes in SpikeBuffer
% Pass them to the saver
% Done !


end

