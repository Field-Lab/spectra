classdef DataFileUpsampler < handle
    %DATAFILEREADER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess = protected, GetAccess = public)
        isBufferUpsampled@logical = false
        isBufferLoaded@logical = false
        isFinished@logical = false
        bufferMaxSize = 4096 % put a power of 2 here for quicker FFTs
        upSampleRatio = 16 % idem
        lastSampleLoaded
        bufferStart
        bufferEnd
        nLPoints
        nRPoints
        
        rawData
        upSampData
        
        nElectrodes
        samplingRate
        
        rawDataFile
        
        startSample
        stopSample
        
        alpha
        
        filterState
        bFilter
        aFilter
    end
    
    methods
        % Constructor
        function obj = DataFileUpsampler(rawDataSource, meanTimeConstant, nLPoints, nRPoints)
            validateattributes(rawDataSource,{'char'},{},'','rawDataSource',1);
            validateattributes(meanTimeConstant,{'numeric'},{'scalar','>',0},'','meanTimeConstant',2);
            validateattributes(nLPoints,{'numeric'},{'scalar','integer','>=',0},'','nLPoints',3);
            validateattributes(nRPoints,{'numeric'},{'scalar','integer','>=',0},'','nRPoints',4);
            
            import edu.ucsc.neurobiology.vision.io.*
            import edu.ucsc.neurobiology.vision.electrodemap.*
            import java.io.*
            
            obj.nLPoints = nLPoints;
            obj.nRPoints = nRPoints;
            
            obj.bufferMaxSize = obj.bufferMaxSize - nLPoints - nRPoints;
            
            parser = DataFileStringParser(rawDataSource);
            datasets = parser.getDatasets();
            obj.rawDataFile = RawDataFile(File(char(datasets(1))));
            startTimes = parser.getStartTimes();
            stopTimes = parser.getStopTimes();
            
            header = obj.rawDataFile.getHeader();
            obj.samplingRate = header.getSamplingFrequency();
            
            obj.startSample = startTimes(1) * obj.samplingRate;
            obj.stopSample = stopTimes(1) * obj.samplingRate;
            obj.lastSampleLoaded = obj.startSample-1;
            
            obj.alpha = 1 / (meanTimeConstant * obj.samplingRate);
            
            packedArrayID = int32(header.getArrayID());
            electrodeMap = ElectrodeMapFactory.getElectrodeMap(packedArrayID);
            obj.nElectrodes = electrodeMap.getNumberOfElectrodes();
            
            obj.filterState = zeros(1,obj.nElectrodes);
            obj.bFilter = (1-obj.alpha)*[1,-1];
            obj.aFilter = [1,obj.alpha-1];
        end
        
        function [bufferStart, bufferEnd] = loadNextBuffer(obj, varargin)
            if nargin == 2 % Can be used to force a different length buffer call.
                validateattributes(varargin{1},{'numeric'},{'scalar','integer','>',0},'','bufferLength',2);
                bufferSize = varargin{1};
            else
                bufferSize = obj.bufferMaxSize;
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
            
            % Filter
            [obj.rawData, obj.filterState] = filter(obj.bFilter, obj.aFilter, ...
                single(obj.rawDataFile.getData(obj.bufferStart, obj.bufferEnd - obj.bufferStart)'),...
                obj.filterState, 2);
            obj.isBufferLoaded = true;
            obj.isBufferUpsampled = false;
            
            % Assign
            bufferStart = obj.bufferStart;
            bufferEnd = obj.bufferEnd;
        end
        
        function upsampleBuffer(obj)
            if ~obj.isBufferLoaded
                throw(MException('','DataFileUpsampler:loadNextBuffer:No Buffer Loaded'));
            end
            if ~obj.isBufferUpsampled
                % rawDataGPU = gpuArray(obj.rawData);
                % fftData = obj.upSampleRatio*fft(rawDataGPU,[],2);
                % obj.upSampData = gather(...
                %   ifft(...
                %   fftData(:,1:(ceil(size(fftData,2)/2))),...
                %   obj.upSampleRatio*size(obj.rawData,2),...
                %   2,'symmetric'));
                fftData = obj.upSampleRatio*fft(obj.rawData,[],2);
                obj.upSampData = ifft(...
                    fftData(:,1:(ceil(size(fftData,2)/2))),...
                    obj.upSampleRatio*size(obj.rawData,2),...
                    2,'symmetric');
            end
            obj.isBufferUpsampled = true;
        end
        
        %         function data = getData(obj,rows,cols)
        
        %         function data = getUpSampData(obj,rows,cols)
    end
    
end

