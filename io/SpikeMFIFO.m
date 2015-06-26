classdef SpikeMFIFO < handle
    %SPIKEMFIFO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable, GetAccess = public)
        timeWindow = 60; % Value after which we are sure no earlier spikes are going to be pushed. Hardcoded 60 samples.
    end
    
    properties (SetAccess = protected, GetAccess = public)
        spikesStored % Currently used spikeArray space
    end
    
    properties (SetAccess = protected, GetAccess = protected)
        spikeArray % Array of spikes sorted by time. A contiguous subpart of the array is used
        firstIndex % First index of spikes currently stored
        lastIndex % Last index of spikes currently stored
        allocSize % Currently allocated spikeArray size
    end
    
    methods
        function obj = SpikeMFIFO()
            obj.spikesStored = 0;
            obj.spikeArray(256,1) = SpikeM(0,0,0); % Initializing array at length 256 with a dummy SpikeM
            obj.firstIndex = 1;
            obj.lastIndex = 0;
            obj.allocSize = 256;
        end % constructor
        
        function val = isEmpty(obj)
            val = obj.spikesStored == 0;
        end
        
        function addLastUnsorted(obj,spikeArray)
            validateattributes(spikeArray,{'SpikeM'},{'nonempty','column'},'','spikeArray');
            n = numel(spikeArray);
            if obj.lastIndex + n < obj.allocSize
                obj.reAlloc(n)
            end
            obj.spikeArray((obj.lastIndex+1):(obj.lastIndex+n)) = spikeArray;
            obj.spikesStored = obj.spikesStored + n;
            obj.lastIndex = obj.lastIndex + n;
        end
        
        function spike = peekFirst(obj)
            spike = obj.spikeArray(obj.firstIndex);
        end
        
        function spike = pollFirst(obj)
            spike = obj.spikeArray(obj.firstIndex);
            obj.firstIndex = obj.firstIndex + 1;
            obj.spikesStored = obj.spikesStored - 1;
        end
        
        function spikes = pollAllOld(obj,nSample) % Remove all spikes older than nSample - timeWindow
            validateattributes(nSample,{'numeric'},{'scalar','integer','>=',0},'','n');
            
            
        end
        
        function reAllocate(obj,n)
            validateattributes(n,{'numeric'},{'scalar','integer','>=',0},'','n');
            % Requires reallocation for insertion of n more objects
            obj.spikeArray(1:obj.spikesStored) = obj.spikeArray(obj.firstIndex:obj.lastIndex);
            obj.firstIndex = 1;
            obj.lastIndex = obj.spikesStored; % Not enough space
            if obj.spikesStored + n > obj.allocSize 
               obj.spikeArray(2^ceil(log2(obj.spikesStored + n))) = SpikeM(0,0,0);
               obj.allocSize = 2^ceil(log2(obj.spikesStored + n));
            end
            if obj.spikesStored + n <= 8 * obj.allocSize % We're having way too much space
               obj.spikeArray = obj.spikeArray(1:obj.spikesStored);
               obj.spikeArray(2^ceil(log2(obj.spikesStored + n))) = SpikeM(0,0,0);
               obj.allocSize = 2^ceil(log2(obj.spikesStored + n));
            end
        end
    end
    
end

