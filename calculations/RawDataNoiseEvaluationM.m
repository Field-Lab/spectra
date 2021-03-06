function rmsNoise = RawDataNoiseEvaluationM(dataPath, saveFolder )
    % RAWDATANOISEEVALUATION Evaluates the amount of noise over all electrodes
    % for a input data recording
    % Generates and save the output .noise file in both ascii single precision
    % number - as column array and .noise.mat file
    %
    % Careful with heading tab on each line generated by matlab save.
    %
    % Inputs:
    %   dataPath: absolute or relative path to data folder/file
    %   saveFolder: absolute or relative path to output save folder
    %
    % Outputs:
    %   rmsNoise: nElectrodes x 1 array containing noise floor values
    %
    % Author -- Vincent Deo -- Stanford University -- August 27, 2015

% Validating attributes
if ~(exist(dataPath,'file') == 2 || exist(dataPath,'file') == 7)
    throw(MException('', 'RawDataNoiseEvaluation: data source does not exist'));
end
if ~(exist(saveFolder,'file') == 7)
    throw(MException('', 'RawDataNoiseEvaluation: dataset folder does not exist'));
end

% Load up config
global GLOBAL_CONFIG
noiseConfig = GLOBAL_CONFIG.getNoiseConfig();

% Set up parameters and data source
time = noiseConfig.time; % Compute noise from 5 secs of data
timeToSkip = noiseConfig.timeToSkip; % Skip the first 5 sec of the file

dataSource = DataFileUpsampler(dataPath);
samplingRate = dataSource.samplingRate;
nSamples = time * samplingRate;

% get data to process noise on
dataSource.loadRandomBuffer(timeToSkip * samplingRate + dataSource.startSample, nSamples, false);


% noise calculation
[rmsNoise,~] = getNoise(dataSource.rawData(2:end,:)); % Remove TTL
rmsNoise = [0;rmsNoise]; % Add back TTL

% Writing output .noise file - as .ascii and as .mat
[~,name,~] = fileparts(dataPath); % Catching folder name as last element of saveFolder path
save([saveFolder,filesep,name,'.noise'],'rmsNoise','-ascii','-double'); % Saving, ascii single precision format
save([saveFolder,filesep,name,'.noise.mat'],'rmsNoise'); % Saving a copy in mat format for further use.

end
