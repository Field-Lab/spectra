classdef DataFileUpsampler < handle
    %DATAFILEREADER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess = public, GetAccess = public)
        % General
        nElectrodes % Number of electrodes
        samplingRate % Sampling rate
        disconnected % disconnected electrodes
        
        % Data buffers
        rawData % nElectrodes x bufferLength buffer
        upSampData % nElectrodes x (bufferLength * upSampRatio) upsampled buffer.
        
        % Cubic Spline interpolation
        interpolant
        isInterpolated@logical = false
        
        % Buffer management
        isBufferUpsampled@logical = false % Upsampled buffer tag. For spike realignment
        isBufferLoaded@logical = false % If the first buffer has ever been loaded
        bufferMaxSize % Maximum buffer size - power of 2 for quick upsampling
        
        lastSampleLoaded % Marker for continous buffers
        bufferStart % Start (inclusive) of current buffer
        bufferEnd % End (exclusive) of current buffer
        
        nLPoints % Number of left point for spike form loading
        nRPoints % Number of right point for spike form loading
        
        % Upsampling management
        upSampleRatio % Precision of upsampling - power of 2 as well
        
        % Data Source management
        rawDataFile % java RawDataFile object
        dataReadThread % java DataFileReadThread object - allows asynchronous ping-pong buffering
        
        startSample % First (inclusive) sample of data source in the data files
        stopSample % Last (exclusive) sample of data source in the data files
        isFinished@logical = false % Tag for last sample of data source reached
        
        % Data Filtering
        alpha % IIR data filter constant
        filterState % Memory of filter state buffer-to-buffer
        bFilter % Filter coefficients
        aFilter % Filter coefficients
    end
    
    methods
        % Constructor
        % Inputs
        % rawDataSource - Data source path + vision style time tags (eg ".../data...(0-10)")
        % (opt) meanTimeConstant - time constant for raw data lowpass filtering
        % (opt) nLPoints - number of points to the left for spike form
        % (opt) nRPoints - number of points to the right for spike form
        %
        % Calls
        % obj = DataFileUpsampler(rawDataSource)
        % obj = DataFileUpsampler(rawDataSource, meanTimeConstant)
        % obj = DataFileUpsampler(rawDataSource, meanTimeConstant, nLPoints, nRPoints)
        function obj = DataFileUpsampler(rawDataSource, varargin)
            narginchk(1,4)
            meanTimeConstant = 1;
            nLPointsInput = 0;
            nRPointsInput = 0;
            
            validateattributes(rawDataSource,{'char'},{},'','rawDataSource',1);
            if nargin >= 2
                meanTimeConstant = varargin{1};
                validateattributes(meanTimeConstant,{'numeric'},{'scalar','>',0},'','meanTimeConstant',2);
            end
            if nargin == 4
                nLPointsInput = varargin{2};
                nRPointsInput = varargin{3};
                validateattributes(nLPointsInput,{'numeric'},{'scalar','integer','>=',0},'','nLPoints',3);
                validateattributes(nRPointsInput,{'numeric'},{'scalar','integer','>=',0},'','nRPoints',4);
            end
            
            import edu.ucsc.neurobiology.vision.io.*
            import edu.ucsc.neurobiology.vision.electrodemap.*
            import java.io.*
            
            obj.nLPoints = nLPointsInput;
            obj.nRPoints = nRPointsInput;
            
            config = mVisionConfig();
            dataConfig = config.getDataConfig();
            obj.bufferMaxSize = dataConfig.bufferMaxSize;
            obj.upSampleRatio = dataConfig.upSampleRatio;
            
            obj.bufferMaxSize = obj.bufferMaxSize - obj.nLPoints - obj.nRPoints;
            
            parser = DataFileStringParser(rawDataSource);
            datasets = parser.getDatasets();
            obj.rawDataFile = RawDataFile(File(char(datasets(1))));
                        
            startTimes = parser.getStartTimes();
            stopTimes = parser.getStopTimes();
            
            header = obj.rawDataFile.getHeader();
            obj.samplingRate = header.getSamplingFrequency();
            
            obj.startSample = startTimes(1) * obj.samplingRate;
            obj.stopSample = min (stopTimes(1) * obj.samplingRate, ...
                header.getNumberOfSamples);
            obj.lastSampleLoaded = obj.startSample-1;
            
            obj.alpha = 1 / (meanTimeConstant * obj.samplingRate);
            
            packedArrayID = int32(header.getArrayID());
            electrodeMap = ElectrodeMapFactory.getElectrodeMap(packedArrayID);
            obj.nElectrodes = electrodeMap.getNumberOfElectrodes();
            obj.disconnected = electrodeMap.getDisconnectedElectrodesList();
           
            obj.dataReadThread = DataFileReadThread(obj.rawDataFile, obj.nElectrodes, obj.alpha);
            obj.dataReadThread.start();
            
            obj.filterState = zeros(1,obj.nElectrodes);
            obj.bFilter = (1-obj.alpha)*[1,-1];
            obj.aFilter = [1,obj.alpha-1];
        end
        
        % Load next buffer in order
        % Maintains the order of the buffers and sequentially loads all the range of the data source
        %
        % Inputs
        % (opt) bufferSize = size of buffer to call. Defaults to the maximal size defined for the class
        % defaults to the size defined as constant class parameter.
        %
        % Calls
        % obj.loadNextBuffer()
        % obj.loadNextBuffer(bufferSize)
        %
        % Returns
        % bufferStart - Inclusive start sample of loaded buffer
        % bufferEnd - Exclusive end sample of loaded buffer
        function [bufferStart, bufferEnd] = loadNextBuffer(obj, varargin)
            narginchk(1,3);
            if nargin >= 2 % Can be used to force a different length buffer call.
                validateattributes(varargin{1},{'numeric'},{'scalar','integer','>',obj.nLPoints+obj.nRPoints},'','bufferLength',2);
                bufferSize = varargin{1}-obj.nLPoints-obj.nRPoints;
                
                % Buffer is not default size. We must clear preloaded default size buffer.
                obj.dataReadThread.clearBuffers();
            else
                bufferSize = obj.bufferMaxSize;
            end
            if nargin == 3
                filterTag = varargin{2};
            else
                filterTag = true; % Default behavior: filter buffers
            end
            
            if obj.isFinished
                throw(MException('','DataFileUpsampler:loadNextBuffer:Reader is already finished'));
            end
            
            obj.bufferStart = max(obj.startSample, obj.lastSampleLoaded + 1 - obj.nLPoints); % Buffer beginning (inclusive)
            obj.bufferEnd = min(obj.lastSampleLoaded + bufferSize + obj.nRPoints + 1, obj.stopSample); % Buffer end (exclusive)
            if obj.bufferEnd == obj.stopSample
                obj.isFinished = true;
            end
            obj.lastSampleLoaded = obj.bufferEnd - obj.nRPoints - 1;
            
            % Collect data
            obj.rawData = obj.dataReadThread.getFilteredBuffer()';
            if numel(obj.rawData) == 0 % First (initialization buffer)
                obj.dataReadThread.loadNextBuffer(obj.bufferStart, obj.bufferEnd - obj.bufferStart, filterTag);
                if filterTag
                    obj.dataReadThread.filterBuffer();
                    obj.rawData = obj.dataReadThread.getFilteredBuffer()';
                else
                    obj.rawData = obj.dataReadThread.getShortBuffer()';
                end               
            end
            
            obj.isBufferLoaded = true;
            obj.isBufferUpsampled = false;
            
            obj.isInterpolated = false;
            
            % Assign
            bufferStart = obj.bufferStart;
            bufferEnd = obj.bufferEnd;
            
            % Compute next buffer parameters and request pong
            nextStart = max(obj.startSample, obj.lastSampleLoaded + 1 - obj.nLPoints); % Buffer beginning (inclusive)
            nextEnd = min(obj.lastSampleLoaded + obj.bufferMaxSize + obj.nRPoints + 1, obj.stopSample); % Buffer end (exclusive)
            % NOTE: ASSUMING NEXT BUFFER WILL BE DEFAULT SIZE AND FILTERED
            
            if ~obj.isFinished
                obj.dataReadThread.loadNextBuffer(nextStart, nextEnd - nextStart, true);
            else
                obj.dataReadThread.setQuit();
            end
            
        end
        
        % Load a requested random access buffer
        % Does not erase the sample markers of buffers loaded with obj.loadNextBuffer()
        % So can be used in the middle of sequential buffer calls
        %
        % Inputs
        % bufferStart = start sample of buffer (inclusive)
        % bufferEnd = stop sample of buffer (exclusive)
        % filterTag = boolean tag if buffer must be filtered.
        % if true, the previous filter state is neither used nor overwritten
        %
        % Calls
        % obj.loadRandomBuffer(bufferStart, bufferEnd, filterTag)
        %
        % Returns
        % bufferStart - Inclusive start sample of loaded buffer
        % bufferEnd - Exclusive end sample of loaded buffer
        function [bufferStart, bufferEnd] = loadRandomBuffer(obj, bufferStart, bufferSize, filterTag)
            validateattributes(bufferStart,{'numeric'},{'scalar','integer','>=',0});
            validateattributes(bufferSize,{'numeric'},{'scalar','integer','>',0},'','bufferLength',2);
            validateattributes(filterTag,{'logical'},{'scalar'},'','filterTag',3);
            
            obj.bufferStart = max(obj.startSample,bufferStart);
            obj.bufferEnd = min(bufferStart + bufferSize, obj.stopSample);
            if obj.bufferEnd == obj.stopSample
                disp('Warning: data source too short to load complete requested random buffer');
            end
            
            % Filter
            if filterTag
                [obj.rawData, ~] = filter(obj.bFilter, obj.aFilter, ...
                    single(obj.rawDataFile.getData(obj.bufferStart, obj.bufferEnd - obj.bufferStart)'),...
                    [], 2);
            else
                obj.rawData = single(obj.rawDataFile.getData(obj.bufferStart, obj.bufferEnd - obj.bufferStart)');
            end
            obj.isBufferLoaded = true;
            obj.isBufferUpsampled = false;
            obj.isInterpolated = false;
            
            % Assign
            bufferStart = obj.bufferStart;
            bufferEnd = obj.bufferEnd;
        end
        
        % Upsampling of currently loaded buffer
        % Method used in non-aliasing constant bandwidth upsampling
        % Instead of cubic spline interpolation as in vision
        %
        % This may be accelerated on GPU
        % Or parallelized by electrode if done in another language
        %
        % Limiting factor is upsampled buffer size in RAM
        % Many other data strategies are worth it
        % this one takes advantage of the large amount of upsampling necessary for all the
        % spikes + Matlab's ability for big ffts.
        % 4096 buffer length * 16 upsampling * 513 electrodes * single precision = 135 MB
        function upsampleBuffer(obj)
            if ~obj.isBufferLoaded
                throw(MException('','DataFileUpsampler:upsampleBuffer:No Buffer Loaded'));
            end
            if ~obj.isBufferUpsampled
                if obj.upSampleRatio > 1
                    fftData = obj.upSampleRatio*fft(obj.rawData,[],2);
                    obj.upSampData = ifft(...
                        fftData(:,1:(ceil(size(fftData,2)/2))),...
                        obj.upSampleRatio*size(obj.rawData,2),...
                        2,'symmetric');
                else
                    obj.upSampData = obj.rawData;
                end
            end
            obj.isBufferUpsampled = true;
        end
        
        function createInterpolant(obj)
            if ~obj.isBufferLoaded
                throw(MException('','DataFileUpsampler:createInterpolant:No Buffer Loaded'));
            end
            
            if ~obj.isInterpolated
                for i = obj.nElectrodes:-1:2;
                    obj.interpolant{i} = griddedInterpolant(1:size(obj.rawData,2),obj.rawData(i,:),'spline');
                end
            end
            
            obj.isInterpolated = true;
        end
        
        % Forces refiltering of current buffer
        % Useful during spikefinding when initial filter state is initialized after analysis of the
        % first second of data
        % Tracking tag sets if the filter state currently stored is to be used and resaved after
        % filtering
        function forceFilter(obj,trackingTag)
            if ~isfloat(obj.rawData)
                obj.rawData = single(obj.rawData);
            end
            
            if trackingTag                
                [obj.rawData, obj.filterState] = filter(obj.bFilter, obj.aFilter, ...
                    obj.rawData, obj.filterState, 2);
                
                obj.dataReadThread.setFilterState(-obj.filterState);
                % We pass the filter state to the IO thread and force refiltering of the NEXT buffer
                % to the one just filtered in Matlab here, which is already preloaded
                obj.dataReadThread.filterBuffer();
            else
                [obj.rawData, ~] = filter(obj.bFilter, obj.aFilter, ...
                    single(obj.rawData),[], 2);
            end
        end
        
    end % methods
end % classdef
