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
            
            validateattributes(timeConstant,{'numeric'},{'scalar','>',0},'','timeConstant');
            obj.alpha = obj.delta / timeConstant;
            
        end % Constructor
        
        
        %%%%%%%%%%%%%%%%% WIP
        function spikes = processBuffer(obj, buffer)
            validateattributes(buffer,{'numeric'},{'2d','nrows',obj.nElectrodes},'','sample');
            if ~obj.initialized
                throw(MException('SpikeFinderM.processSample: not initialized'));
            end
            
            % processing of TTL buffer
            ttlThresholded = [obj.buildingSpike(1),buffer(1,:) < -obj.ttlThreshold];
            ttlThUp = find(and(~ttlThresholded(1:(end-1)),ttlThresholded(2:end)));
            ttlThDown = find(and(ttlThresholded(1:(end-1)),~ttlThresholded(2:end)));
            
            ttlSpike = [ttlThUp',repmat([1,1500],numel(ttlThUp),1)];
            
            if numel(ttlThUp) > 0
                if obj.previousSpikeTime <= 0 % First ttl ever is in this buffer
                    obj.ttlAverage = obj.ttlAverage + ttlThUp(end) - ttlThUp(1);
                    obj.ttlCount = obj.ttlCount + numel(ttlThUp) - 1;
                else
                    obj.ttlAverage = obj.ttlAverage + ttlThUp(end) - obj.previousSpikeTime(1);
                    obj.ttlCount = obj.ttlCount + numel(ttlThUp);
                end
                
                obj.previousSpikeTime(1) = ttlThUp(end);
                
                if numel(ttlThDown) > 0
                    obj.buildingSpike(1) = ttlThDown(end) < ttlThUp(end);
                else
                    obj.buildingSpike(1) = true;
                end
            else
                obj.buildingSpike(1) = obj.buildingSpike(1) && numel(ttlThDown) == 0;
            end
            
            spikes = ttlSpike;
            %%%
            % processing of all other electrodes
            bufferThresholded = [obj.buildingSpike(1:end),bsxfun(@lt,buffer(1:end,:),-obj.spikeThresholds)];
            
            for el = 2:obj.nElectrodes
                if obj.disconnected(el)
                    continue
                end
                
                elThUp = find(and(~bufferThresholded(el,1:(end-1)),bufferThresholded(el,2:end)));
                elThDown = find(and(bufferThresholded(el,1:(end-1)),~bufferThresholded(el,2:end)));
                
                % Align frames
                if numel(elThUp) > numel(elThDown)
                    frames = [elThUp(1:(end-1)) ; elThDown];
                    lastUp = elThUp(end);
                    firstDown = [];
                else if numel(elThUp) < numel(elThDown)
                        frames = [elThUp(1:end) ; elThDown(2:end)];
                        lastUp = [];
                        firstDown = elThDown(1);
                    else
                        if numel(elThUp) == 0
                            continue
                        end
                        if elThUp(1) < elThDown(1)
                            frames = [elThUp ; elThDown];
                            lastUp = [];
                            firstDown = [];
                        else
                            frames = [elThUp(1:(end-1));elThDown(2:end)];
                            firstDown = elThDown(1);
                            lastUp = elThUp(end);
                        end
                    end
                end
                
                spikesEl = zeros(0,3);
                % check the case where first, last and frames ar empty
                
                if numel(firstDown) > 0
                    x = firstDown;
                    
                    time = find(buffer(el,1:x) == min(buffer(el,1:x)),1);
                    amp = buffer(el,time);
                    if amp < obj.maxAmplitude(el);
                        obj.maxTime(el) = time + obj.currentSample;
                        obj.maxAmplitude(el) = amp;
                    end
                    
                    if (x + obj.currentSample - obj.startTime(el) <= obj.maxSpikeWidth) &&...
                            (obj.maxTime(el) - obj.previousSpikeTime(el) > obj.minTimeSeparation)
                        % Spike is valid inter-time-wise
                        spikesEl = [spikesEl;obj.maxTime(el),el,obj.maxAmplitude(el)];
                    end
                    
                    obj.buildingSpike(el) = false;
                    obj.previousSpikeTime(el) = obj.maxTime(el);
                    
                else if numel(frames) == 0 && numel(lastUp) == 0 && obj.buildingSpike(el)
                        % Corner case - spikes running across whole buffer
                        time = find(buffer(el,:) == min(buffer(el,:)),1);
                        amp = buffer(el,time);
                        if amp < obj.maxAmplitude(el);
                            obj.maxTime(el) = time + obj.currentSample;
                            obj.maxAmplitude(el) = amp;
                        end
                    end
                end
                
                for f = frames
                    time = f(1) - 1 + find(buffer(el,f(1):f(2)) == min(buffer(el,f(1):f(2))),1);
                    amp = buffer(el,time);
                    
                    if ((f(2) - f(1)) <= obj.maxSpikeWidth &&...
                       (time + obj.currentSample - obj.previousSpikeTime(el)) > obj.minTimeSeparation)
                        % Spike is valid inter-time-wise
                        spikesEl = [spikesEl;time + obj.currentSample,el,amp];
                    end
                    
                    obj.previousSpikeTime(el) = time + obj.currentSample;
                    
                end
                
                if numel(lastUp) > 0
                    x = lastUp;
                    obj.maxTime(el) = find(buffer(el,x:end) == min(buffer(el,x:end)),1) + obj.currentSample;
                    obj.maxAmplitude(el) = buffer(el,time);
                    
                    obj.buildingSpike(el) = true;
                end
                
                spikes = [spikes;spikesEl];
            end
            
            spikes = sortrows(spikes,1);
            
            % updating spikes counter
            obj.totalSpikes = obj.totalSpikes + size(spikes,1);
            
            % updating current sample
            obj.currentSample = obj.currentSample + size(buffer,2);
        end
        %%%%%%%%%%%%%%%%%%%
        
        function spikes = processSample(obj, sample)
            validateattributes(sample,{'numeric'},{'column','nrows',obj.nElectrodes},'','sample');
            if ~obj.initialized
                throw(MException('SpikeFinderM.processSample: not initialized'));
            end
            
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
            
            % Filtering is done in dataFileUpsampler - we just flip sign here
            amplitude = -sample;
            
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
            obj.totalSpikes = obj.totalSpikes + size(spikes,1);
            
            % updating current sample
            obj.currentSample = obj.currentSample + 1;
            
        end % processSample
        
        function meanCorrection = initialize(obj,varargin)
            narginchk(1,2);
            if obj.initialized
                throw(MException('SpikeFinderM_initialize:InitializationError','Multiple Initialization'));
            end
            if obj.currentSample ~= 0
                throw(MException('SpikeFinderM_initialize:InitializationError','Sample counter already > 0'));
            end
            if nargin == 1
                obj.initialize(zeros(obj.nElectrodes,0));
                obj.initialized = true;
            else
                sampleBuffer = varargin{1};
                validateattributes(sampleBuffer,{'numeric'},{'2d','nrows',obj.nElectrodes},'','sampleBuffer');
                meanCorrection = mean(sampleBuffer,2);
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
    
end