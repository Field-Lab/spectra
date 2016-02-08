% Grab -any- header ! (just with the right number of electrodes necessarily, arrayID optionally)
% Note: check your java path
%%%%%

% Make a header
secondsTime = 123456789;
nElectrodes = 2;
frequency = 20000;
nSamples = 233600;
arrayID = 9999;
format = 1;
datasetIdentifier = 'data999-data999';
comment = '';

hdr = edu.ucsc.neurobiology.vision.io.RawDataHeader512(...
    secondsTime, nElectrodes, frequency,...
    nSamples, arrayID,...
    format, datasetIdentifier, comment);

% no need to touch nSamples, that's automatic at read.
% ---> THIS HAS TO BE data<3digits>-data<same3digits>. Because vision.

% Instantiate a raw data saver.
saver = edu.ucsc.neurobiology.vision.io.RawDataSaver(...
    '2005-04-26-0',... % don't leave empty string (java gets a null because matlab is a jackass)
    'X:\EJgroup_data\data',... % here neither
    hdr,...
    2048,...
    4,...
    secondsTime);

%%
rdf = edu.ucsc.neurobiology.vision.io.RawDataFile('X:\EJgroup_data\data\2005-04-26-0\data002');
myCustomData = rdf.getData(0,233610);
myCustomData = myCustomData(:,1:2);
%%
for i = 1:nSamples %nSamples
    if mod(i,20000) == 0
        fprintf('%u\n',i);
    end
    saver.processSample(myCustomData(i,:));
end
saver.finishSampleProcessing();
