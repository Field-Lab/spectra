classdef SpikeFinderM < handle
    %SPIKEFINDER Encapsulating the process of finding, creating and sorting spikes
    %   Replaces the spikeFinder + spikeBuffer in the java architecture
    %   On the concept of listener called -> call listeners
    %   But processed on sample arrays.
    properties (SetAccess = immutable, GetAccess = protected)
        delta = 50e-6     % Sampling timestep
        minTimeSeparation % .25 msec separation, in samples
        maxSpikeWidth     % 2.5 msec maximum spike width, in samples
        
        nElectrodes
        disconnected @logical
        spikeThresholds
        ttlThreshold
        alpha
        dataFileUpsampler
    end
    properties (SetAccess = protected, GetAccess = public)
        buildingSpike @logical
        initialized = false
        maxTime
        maxAmplitude
        previousSpikeTime
        startTime
        currentSample = 0; % Last sample of the last treated buffer % Unused if per-buffer finding
        ttlAverage = 0;
        ttlCount = 0; % Counter for ttl intervals. 1 less than TTL spikes outputed
        totalSpikes = 0;
    end
    
    methods
        % Constructor
        function obj = SpikeFinderM(spikeThresholds, ttlThreshold, timeConstant, dataFileUpsampler)
            
            validateattributes(ttlThreshold,{'numeric'},{'scalar'},'','ttlThresholds');
            obj.ttlThreshold = double(ttlThreshold);
            
            validateattributes(dataFileUpsampler,{'DataFileUpsampler'},{},'','dataFileUpsampler',5);
            obj.dataFileUpsampler = dataFileUpsampler;
            obj.nElectrodes = dataFileUpsampler.nElectrodes;
            obj.disconnected = dataFileUpsampler.disconnected;
            
            validateattributes(spikeThresholds,{'numeric'},{'column','nrows',obj.nElectrodes},'','spikeThresholds');
            obj.spikeThresholds = double(spikeThresholds);
            
            % For process flow
            obj.buildingSpike = false(obj.nElectrodes, 1);
            obj.maxTime = zeros(obj.nElectrodes,1);
            obj.previousSpikeTime = -1000*ones(obj.nElectrodes,1);
            obj.maxAmplitude = zeros(obj.nElectrodes,1);
            obj.startTime = zeros(obj.nElectrodes,1);
            
            config = mVisionConfig();
            spConfig = config.getSpikeConfig();
            obj.maxSpikeWidth = spConfig.maxSpikeWidth;
            obj.minTimeSeparation = spConfig.minSpikeSeparation;
            
            validateattributes(timeConstant,{'numeric'},{'scalar','>',0},'','timeConstant');
            obj.alpha = obj.delta / timeConstant;
            
        end % Constructor
        
        
        % Processing of a buffer of data nElectrodes * nSamples
        % And outputs all the finished spikes, in order
        % No use of spike buffer anymore: spikes on buffer overlap may not be sorted
        % Final sort required before saving.
        function spikes = processBuffer(obj)
            
            bs = obj.dataFileUpsampler.bufferStart;
            
            spikeStore = cell(obj.nElectrodes,1);
            
            if ~obj.initialized
                throw(MException('SpikeFinderM.processSample: not initialized'));
            end
            
            % processing of TTL buffer
            ttlThresholded = [obj.buildingSpike(1),obj.dataFileUpsampler.rawData(1,:) < -obj.ttlThreshold];
            ttlThUp = find(and(~ttlThresholded(1:(end-1)),ttlThresholded(2:end)));
            ttlThDown = find(and(ttlThresholded(1:(end-1)),~ttlThresholded(2:end)));
            
            ttlSpike = [ttlThUp'+bs-1,repmat([1,1500],numel(ttlThUp),1)];
            
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
            
            spikeStore{1} = ttlSpike;
            
            %%%
            % processing of all other electrodes
            bufferThresholded = [obj.buildingSpike(1:end),bsxfun(@lt,obj.dataFileUpsampler.rawData(1:end,:),-obj.spikeThresholds)];
            
            xorBuff = xor(bufferThresholded(:,2:end),bufferThresholded(:,1:(end-1)));
            
            % Allocation - extension of num rows done further down if needed
            frameStack = nan(1,obj.maxSpikeWidth + 1);
            
            for el = (find(~obj.disconnected(2:end))+1)'
                
                elThCross = find(xorBuff(el,:));
                
                try
                    if bufferThresholded(el,elThCross(1)+1);
                       elThUp = elThCross(1:2:end);
                       elThDown = elThCross(2:2:end);
                    else
                       elThUp = elThCross(2:2:end);
                       elThDown = elThCross(1:2:end);
                    end
                catch e
                    elThUp = [];
                    elThDown = [];
                end
                
                
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
                    
                    [amp,time] = min(obj.dataFileUpsampler.rawData(el,1:x));
                    
                    if amp < obj.maxAmplitude(el);
                        obj.maxTime(el) = time + bs - 1;
                        obj.maxAmplitude(el) = amp;
                    end
                    
                    if (x + bs - obj.startTime(el) <= obj.maxSpikeWidth) &&...
                            (obj.maxTime(el) - obj.previousSpikeTime(el) > obj.minTimeSeparation)
                        % Spike is valid inter-time-wise
                        spikesEl = [spikesEl;obj.maxTime(el),el,-obj.maxAmplitude(el)];
                    end
                    
                    obj.buildingSpike(el) = false;
                    obj.previousSpikeTime(el) = obj.maxTime(el);
                    
                else if numel(frames) == 0 && numel(lastUp) == 0 && obj.buildingSpike(el)
                        % Corner case - spikes running across whole buffer
                        [amp,time] = min(obj.dataFileUpsampler.rawData(el,:));
                        if amp < obj.maxAmplitude(el);
                            obj.maxTime(el) = time + bs - 1;
                            obj.maxAmplitude(el) = amp;
                        end
                    end
                end
                
                if numel(frames) > 0
                    
                    frameLength = frames(2,:) - frames(1,:);
                    
                    if size(frames,2) > size(frameStack,1) || (max(frameLength)+1) > size(frames,1) 
                        frameStack = Inf(max(size(frames,2),size(frameStack,1)),max(max(frameLength)+1,size(frames,1)));
                    else
                        frameStack(:) = nan;
                    end
                    
                    buffSeed = zeros(1,sum(frameLength)+ size(frameLength,2));
                    frameDest = zeros(1,sum(frameLength)+ size(frameLength,2));
                    ind = [1,cumsum(frameLength+1)+1];
                    
                    for f = 1:size(frames,2)
                        buffSeed(ind(f):(ind(f)+frameLength(f))) = (0:frameLength(f)) + frames(1,f);
                        frameDest(ind(f):(ind(f)+frameLength(f))) = (0:frameLength(f)) * size(frameStack,1) + f;
                    end
                    
                    frameStack(frameDest) = obj.dataFileUpsampler.rawData(el,buffSeed);
               
                    [amp,I] = min(frameStack,[],2,'omitnan');
                    time = frames(1,:) - 1 + I(1:size(frames,2))';
                    
                    timeSpacing = time - [obj.previousSpikeTime(el) - bs + 1, time(1:(end-1))];
                    keep = and(frameLength <= obj.maxSpikeWidth, timeSpacing > obj.minTimeSeparation);

                    spikesEl = [spikesEl; [time(keep)'+ bs - 1,repmat(el,nnz(keep),1),-amp(keep)]];
            
                end
                
                if numel(lastUp) > 0
                    x = lastUp;
                    [amp,time] = min(obj.dataFileUpsampler.rawData(el,x:end));
                    obj.maxTime(el) = time + bs + x - 2;
                    obj.maxAmplitude(el) = amp;
                    
                    obj.startTime(el) = x + bs - 1;
                    
                    obj.buildingSpike(el) = true;
                end
                
                spikeStore{el} = spikesEl;
            
                if false % Alignment debug plots
                    %%
                    clf
                    plot(obj.dataFileUpsampler.rawData(el,:),'b+-');
                    hold on
                    plot([0,size(obj.dataFileUpsampler.rawData,2)],-obj.spikeThresholds(el)*[1,1],'r-');
                    stem(spikesEl(:,1)-bs,-spikesEl(:,3),'r');
                    hold off
                end
            
            end
            
            % Merging and sorting spikes
            spikes = vertcat(spikeStore{:});
            spikes = sortrows(spikes,1);
            
            % Updating spikes counter
            obj.totalSpikes = obj.totalSpikes + size(spikes,1);
            
        end % processBuffer
        
        
        % ---- OBSOLETE ----
        % Processing of a single sample over nElectrodes
        % And outputs the spikes terminating at this sample.
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
            if nargin == 1 % Initialize on datasource buffer
                meanCorrection = mean(obj.dataFileUpsampler.rawData,2);
                obj.initialized = true;
            else % initialize on argument buffer
                validateattributes(varargin{1},{'numeric'},{'2d','nrows',obj.nElectrodes},'','sampleBuffer');
                meanCorrection = mean(varargin{1},2);
                obj.initialized = true;
            end % nargin test
        end % initialize
        
    end
    
end