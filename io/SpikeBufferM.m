classdef SpikeBufferM < handle
    %SPIKEBUFFERM covers storing and sorting of spikes found by the spikeFinder
    % Stores a SpikeM array, keeps it sorted by time
    
    properties (SetAccess = immutable, GetAccess = public)
        timeWindow = 60; % Value after which we are sure no earlier spikes are going to be pushed. Hardcoded 60 samples.
    end
    
    properties (SetAccess = protected, GetAccess = protected)
        spikeArray % Array of spikes sorted by time. A contiguous subpart of the array is used
    end
    
    methods
        function obj = SpikeBufferM()
           obj.spikeArray = []; 
        end
        
        function addSpikes(obj,spikes)
           validateattributes(spikes,{'numeric'},{'2d','ncols',3},'','spikes');
           obj.spikeArray = [obj.spikeArray; spikes];
        end
        
        function spikes = getSpikes(obj,nSample)
            validateattributes(nSample,{'numeric'},{'scalar'},'','nSample');
            spikes = obj.spikeArray(nSample - obj.spikeArray(:,1) >= obj.timeWindow, :);
            spikes = sortrows(spikes,1);
            obj.spikeArray = obj.spikeArray(~(nSample - obj.spikeArray(:,1) >= obj.timeWindow), :);
        end
        
        function spikes = getAllSpikes(obj)
            spikes = sortrows(obj.spikeArray,1);
            obj.spikeArray = [];
        end
    end
    
end

