classdef SpikeFinderM < handle
    %SPIKEFINDER Encapsulating the process of finding, creating and sorting spikes
    %   Replaces the spikeFinder + spikeBuffer in the java architecture
    %   On the concept of listener called -> call listeners
    %   But processed on sample arrays.
    properties (Constant, Access = protected)
        delta = 50e-6; % Sampling timestep
        minTimeSeparation = 0.25 * 20; % .25 msec separation, in samples
        maxSpikeWidth = 2.5 * 20; % 2.5 msec maximum spike width, in samples
    end
    properties (SetAccess = immutable, GetAccess = public)
        nElectrodes
        disconnected @logical
        spikeThresholds
        ttlThreshold
        alpha
    end
    properties (SetAccess = protected, GetAccess = public)
        buildingSpike @logical
        initialized = false
        maxTime
        maxAmplitude
        previousSpikeTime
        startTime
        filterState
        currentSample = 0; % Last sample of the last treated buffer
        ttlAverage = 0;
        ttlCount = 0; % Counter for ttl intervals. 1 less than TTL spikes outputed
        totalSpikes = 0;
    end
    
    methods
        % Constructor
        function obj = SpikeFinderM(electrodeMap, spikeThresholds, ttlThreshold, timeConstant)
            % validate electrode map
            validateattributes(electrodeMap,{'edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMap'},{},'','electrodeMap');
            obj.nElectrodes = int32(electrodeMap.getNumberOfElectrodes());
            obj.disconnected = electrodeMap.getDisconnectedElectrodesList();
            
            validateattributes(spikeThresholds,{'numeric'},{'column','nrows',obj.nElectrodes},'','spikeThresholds');
            obj.spikeThresholds = double(spikeThresholds);
            
            validateattributes(ttlThreshold,{'numeric'},{'scalar'},'','ttlThresholds');
            obj.ttlThreshold = double(ttlThreshold);
            
            % For process flow
            obj.buildingSpike = false(obj.nElectrodes, 1);
            obj.maxTime = zeros(obj.nElectrodes,1);
            obj.previousSpikeTime = -1000*ones(obj.nElectrodes,1);
            obj.maxAmplitude = zeros(obj.nElectrodes,1);
            obj.startTime = zeros(obj.nElectrodes,1);
            obj.filterState = zeros(obj.nElectrodes,1);
            
            validateattributes(timeConstant,{'numeric'},{'scalar','>',0},'','timeConstant');
            obj.alpha = obj.delta / timeConstant;
            
        end % Constructor
        
        function spikes = processSample(obj, sample)
            validateattributes(sample,{'numeric'},{'column','nrows',obj.nElectrodes},'','sample');
            if ~obj.initialized
                throw(MException('SpikeFinderM.processSample: not initialized'));
            end
            % updating current sample
            obj.currentSample = obj.currentSample + 1;
            
            % processing of TTL sample
            ttlSpike = zeros(0,3);
            if sample(1) < -obj.ttlThreshold
                if ~obj.buildingSpike(1)
                    obj.buildingSpike(1) = true;
                    ttlSpike = [obj.currentSample,1,1500]; %% Hardcoded 1500 amplitude for ttl spikes
                    if obj.previousSpikeTime(1) > 0
                        obj.ttlAverage = obj.ttlAverage + obj.currentSample - obj.previousSpikeTime(1);
                        obj.ttlCount = obj.ttlCount + 1;
                    end
                    obj.previousSpikeTime(1) = obj.currentSample;
                end
            else
                obj.buildingSpike(1) = false;
            end
            
            % processing of all other electrodes
            
            % filtering
            obj.filterState(2:end) = (1-obj.alpha) * obj.filterState(2:end) + obj.alpha * sample(2:end);
            amplitude = obj.filterState - sample;
            % % Non filtered
            % amplitude = -sample;
            
            % Identifying behaviors and removing disconnected electrodes
            % Spikes starting at this sample
            startSpike = and(~obj.disconnected,...
                and(amplitude > obj.spikeThresholds, ~obj.buildingSpike));
            % Spikes already started still above threshold
            contSpike = and(~obj.disconnected,...
                and(and(amplitude > obj.spikeThresholds, obj.buildingSpike),...
                amplitude > obj.maxAmplitude));
            % Spikes that just went below samples
            terminateSpike = and(~obj.disconnected,...
                and(amplitude <= obj.spikeThresholds, obj.buildingSpike));
            % Spikes that just went below sample AND are valide - we create a spike from those.
            closeSpike = and(obj.currentSample - obj.startTime <= obj.maxSpikeWidth,...
                and(obj.maxTime - obj.previousSpikeTime > obj.minTimeSeparation,...
                terminateSpike));
            
            % Eliminating TTL electrode from the above filters
            startSpike(1) = false;
            contSpike(1) = false;
            terminateSpike(1) = false;
            closeSpike(1) = false;
            
            % Starting spikes
            if any(startSpike)
                obj.maxTime(startSpike) = obj.currentSample;
                obj.startTime(startSpike) = obj.currentSample;
                obj.maxAmplitude(startSpike) = amplitude(startSpike);
                obj.buildingSpike(startSpike) = true;
            end
            
            % Continuing and rising spikes
            if any(contSpike)
                obj.maxTime(contSpike) = obj.currentSample;
                obj.maxAmplitude(contSpike) = amplitude(contSpike);
            end
            
            % Terminating spikes and assigning output
            spikes = ttlSpike;
            
            if any(terminateSpike)
                obj.previousSpikeTime(terminateSpike) = obj.maxTime(terminateSpike);
                obj.buildingSpike(terminateSpike) = false;
                
                if any(closeSpike)
                    electrode = (1:obj.nElectrodes)';
                    spikes = [ttlSpike;[obj.maxTime(closeSpike),electrode(closeSpike),obj.maxAmplitude(closeSpike)]];
                end
            end
            %             spikes
            %             numel(spikes)
            obj.totalSpikes = obj.totalSpikes + numel(spikes);
            
        end % processSample
        
        function initialize(obj,varargin)
            narginchk(1,2);
            if obj.initialized
                throw(MException('SpikeFinderM_initialize:InitializationError','Multiple Initialization'));
            end
            if obj.currentSample ~= 0
                throw(MException('SpikeFinderM_initialize:InitializationError','Sample counter already > 0'));
            end
            if nargin == 1
                obj.initialize(zeros(513,0));
                obj.initialized = true;
            else
                sampleBuffer = varargin{1};
                validateattributes(sampleBuffer,{'numeric'},{'2d','nrows',obj.nElectrodes},'','sampleBuffer');
                obj.filterState = mean(sampleBuffer,2);
                obj.initialized = true;
            end % nargin test
        end % initialize
        
        function rate = getRefreshRate(obj)
            rate = obj.ttlAverage/obj.ttlCount;
        end % getRefreshRate
        
        function samples = getSamplesProcessed(obj)
            samples = obj.currentSample;
        end % getSamplesProcessed
        
        function spikes = getSpikesFound(obj)
            spikes = obj.totalSpikes;
        end % getSpikesFound
        
    end
    
    % % Draft tentative methods - goal is to process multiple samples at the same time over
    % % several/multiple electrodes.
    %     methods %(Access = protected)
    %         function spikeArray = processTtlBuffer(obj, sampleBuffer)
    %             % Detecting TTL thresholds
    %             if ~obj.buildingSpike(1) && sampleBuffer(1) < -obj.ttlThreshold % first sample of buffer is a new spike
    %                 upFronts = [1,strfind(sampleBuffer < -obj.ttlThreshold,[0 1])+1];
    %             else
    %                 upFronts = strfind(sampleBuffer < -obj.ttlThreshold,[0 1])+1; % strfind yields index of [0 1] pattern, we want the index of the 1.
    %             end
    %
    %             % Creating TTL spikes
    %             nSpikes = size(upFronts,2);
    %             for n = nSpikes:-1:1
    %                 spikeArray(1,n) = SpikeM(obj.currentSample+upFronts(n),1,1500);
    %             end
    %
    %             % Updating finder parameters for next buffer
    %             obj.buildingSpike(1) = sampleBuffer(end) < -obj.ttlThreshold;
    %             if obj.previousSpikeTime(1) > 0
    %                 obj.ttlAverage = obj.ttlAverage + spikeArray(1,n).time - obj.previousSpikeTime(1);
    %                 obj.ttlCount = obj.ttlCount + nSpikes;
    %             end
    %             obj.previousSpikeTime = spikeArray(1,n).time;
    %             obj.totalSpikes = obj.totalSpikes + nSpikes;
    %
    %         end % processTtlBuffer
    %
    %         function spikeArray = processElectrodeBuffer(obj, el, sampleBuffer)
    %             % Filter data
    %             dataFilt = filter(obj.alpha,[1,obj.alpha-1],sampleBuffer,obj.filterState(el));
    %             % Memorizing new filter state
    %             obj.filterState(el) = dataFilt(end);
    %             % Substracting filtered central value - invert sign
    %             data = dataFilt - sampleBuffer;
    %
    %             % Finding up and down fronts along the buffer
    %             % adding an up front at position 1 if a spike appears at this very sample
    %             % upFront values correspond to the first 1 of the spike
    %             % downFront values correspond to the first 0 after the spike
    %             aboveThreshold = data > obj.spikeThresholds(el);
    %             if ~obj.buildingSpike(el) && aboveThreshold(1)
    %                 upFronts = [1,strfind(aboveThreshold,[0 1])+1];
    %             else
    %                 upFronts = strfind(aboveThreshold,[0 1])+1;
    %             end
    %             if obj.buildingSpike(el) && ~aboveThreshold(1)
    %                 downFronts = [1,strfind(aboveThreshold,[1 0])+1];
    %             else
    %                 downFronts = strfind(aboveThreshold,[1 0])+1;
    %             end
    %
    %             % Creating all finished spikes EXCEPT the first one
    %             nSpikes = size(downFronts,2);
    %             for n = nSpikes:-1:2
    %                 maxSpikeTime(1,n) = 0;
    %
    %                 spikeArray(1,n) = SpikeM(obj.currentSample+maxSpikeTime(n),el,maxAmplitude);
    %             end
    %
    %
    %             if obj.buildingSpike(el) % Case: we need to finish building a running spike
    %                 if isempty(downFronts) % No downfronts - spike runs across the whole buffer
    %
    %                 else % Current spike runs up to the first down front
    %
    %                 end
    %             end
    %         end % processElectrodeBuffer
    
    %     end
end