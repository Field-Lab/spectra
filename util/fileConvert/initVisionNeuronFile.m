function newNeuronFile = initVisionNeuronFile(rawDataFilePath, newNeuronFilePath, ttlTimes, varargin)
%INITVISIONNEURONFILE Initializes a Vision neuron file.
% 
%   INITVISIONNEURONFILE(RAWDATAFILEPATH, NEWNEURONFILEPATH, TTLTIMES, ...
%      VISIONPATH) 
%   abstracts away some of the gritty details associated with
%   creating a neuron file that makes sense to Vision.
%   RAWDATAFILEPATH should be the path to a Vision bin file with a header
%   or a folder of bin files. NEWNEURONFILEPATH is the path where the new
%   neuron file will be created. TTLTIMES is a vector of TTL times as 
%   recorded on channel 0 of the raw data file. 
%   VISIONPATH is optional and can be specified if the Vision jar file 
%   has not yet been added to the java class path.
%
%   The function returns a edu.ucsc.neurobiology.vision.io.NeuronFile object
%   to which neurons can be added using the class method:
%      neuronFile.addNeuron(int electrode, int neuronID, int[] times, 
%          int nSpikes)
%   Note the integer type of the arguments.
%   If you don't know what your neuronID is supposed to be, you can call
%      neuronFile.getNeuronID(electrode, clusterIndex) 
%   and Vision will calculate it for you.
%
%   Don't forget to close() the neuronFile after you're done with it.
%   

% ------------
% Add - Vincent Deo - Stanford University - 06/02/2015
% If the jar is already imported in call hierarchy, any argument can be passed as
% visionPath. Changed 4th argument to optional.
narginchk(3,4);
if nargin == 4
    visionPath = varargin{1};
end
% -------------

% Check we linked to the Vision jar file already
if ~exist('edu/ucsc/neurobiology/vision/io/RawDataFile','class') && exist('visionPath', 'var')
    javaaddpath(visionPath);
end

% Get the bin file header 
rawDataFile = edu.ucsc.neurobiology.vision.io.RawDataFile(rawDataFilePath);
rawDataHeader = rawDataFile.getHeader();
rawDataFile.close();

% Magic number for the Vision Header. 
% If things break look for the value of 
%   edu.ucsc.neurobiology.vision.util.VisionParams.NEURONS_HEADER_CAPACITY
% in the Vision source code. 
NEURONS_HEADER_CAPACITY = edu.ucsc.neurobiology.vision.util.VisionParams.NEURONS_HEADER_CAPACITY;
% Same with the Vision Version value (also found in VisionParams)
VISION_VERSION = edu.ucsc.neurobiology.vision.util.VisionParams.VERSION;
% Now the neurons file header version is hidden in:
%   edu.ucsc.neurobiology.vision.io.NeuronFile.INT_VERSION
VERSION = edu.ucsc.neurobiology.vision.io.NeuronFile.INT_VERSION;

% Instantiate a new Vision header 
visionHeader = edu.ucsc.neurobiology.vision.io.VisionHeader();

% There's of course a MAGIC number because Vision likes magic.
% As far as I can tell it's supposed to be:
%   edu.ucsc.neurobiology.vision.io.ProjectionFile.MAGIC
% which is currently set to 0xBBBBBB
MAGIC = hex2dec('BBBBBB');

% Then some reasonable defaults. Not sure why they're needed but 
% Vision specifies them and it's not a good idea to not specify what 
% Vision specifies, in general, because then things break in unexpected
% places.
MIN_SPIKES = 100;
MAX_CONTAM = 1;

% Fill in the Vision header object
visionHeader.magic = MAGIC;
visionHeader.headerVersion = 1;
visionHeader.version = VERSION;
visionHeader.meanTimeConstant = -1;
visionHeader.threshold = -1;
visionHeader.arrayID = rawDataHeader.getArrayID();
visionHeader.nSamples = rawDataHeader.getNumberOfSamples();
visionHeader.samplingFrequency = rawDataHeader.getSamplingFrequency();
visionHeader.visionVersion = VISION_VERSION;
visionHeader.nDimensions = 5;
visionHeader.minNeuronSpikes = MIN_SPIKES;
visionHeader.maxContamination = MAX_CONTAM;

% TTL times should be integers. Cast them to int32 just in case
% Vision neuron files also don't like empty TTL lists, so if there were none,
% we just add a TTL at time 0
if isempty(ttlTimes)
    ttlTimes = int32([0]);
else
    ttlTimes = int32(ttlTimes);
end

% Now we can instantiate the new neuron file
newNeuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(newNeuronFilePath, ...
    visionHeader, NEURONS_HEADER_CAPACITY, ttlTimes);

end % initVisionNeuronFile
