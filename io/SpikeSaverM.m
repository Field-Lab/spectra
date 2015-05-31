classdef SpikeSaverM
    %SPIKESAVERM Wrapper of the java class SpikeSaver
    % Wraps a private SpikeSaver (java) object.
    % Callers to useful java methods are implemented a Matlab methods
    % Along with required getters
    
    properties (SetAccess = immutable, GetAccess = protected)
        % Put some java object here
        % TODO change GetAccess to protected and implement some proper getters and setter
        
        spikeSaver % @edu.ucsc.neurobiology.vision.io.SpikeSaver - built in constructor
        spikeFile % @edu.ucsc.neurobiology.vision.io.SpikeFile - built in constructor
    end
    properties (SetAccess = immutable, GetAccess = public)
        savePath@char % Output folder in String format - provided to constructor and checked
    end
    
    methods
        % put getters, spike treater & debug stuff
        function obj = SpikeSaverM(header, savePath, meanTimeConstant, spikeThreshold)
            import edu.ucsc.neurobiology.vision.io.*
            
            % header class check validity check
            validateattributes(header, {'edu.ucsc.neurobiology.vision.io.RawDataHeader512'},{},'','header',1);
            
            % TODO Complete a path validity check for outputPath argument
            % savePath validity check
            validateattributes(savePath,{'char'},{},'','savePath',2);
            obj.savePath = savePath;
            
            % Create spikeFile object. outputPath should be corrected for
            % the argument should be a folder path and is used here as a file path
            
            % header does NOT always contain the correct packed array ID. See SpikeFinding.java
            
            [~,name,~] = fileparts(savePath);
            obj.spikeFile = SpikeFile([savePath,filesep,name,'.spikes'], header.getArrayID(),...
                meanTimeConstant, spikeThreshold,...
                header.getNumberOfSamples(),...
                header.getSamplingFrequency());
            
            obj.spikeSaver = SpikeSaver(obj.spikeFile);
            
        end
        
        function spikesFound = getSpikesFound(obj)
            spikesFound = obj.spikeSaver.getSpikesFound();
        end
        
        function processSpikes(obj, spikes)
            import edu.ucsc.neurobiology.vision.io.*
            
            validateattributes(spikes,{'numeric'},{'2d','ncols',3},'','spikes');
            validateattributes(spikes(:,1),{'numeric'},{'nondecreasing','integer'},'','spikes(:,1)');
            
            for i = 1:size(spikes,1)
                spikeJava = Spike(spikes(i,1),spikes(i,2)-1,spikes(i,3));
                obj.spikeSaver.processSpike(spikeJava);
            end
        end
        
        function finishSpikeProcessing(obj)
            obj.spikeSaver.finishSpikeProcessing();
        end
    end
    
end

