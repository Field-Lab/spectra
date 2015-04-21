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
        function obj = SpikeSaverM(sampleInputStream, savePath)
            import edu.ucsc.neurobiology.vision.io.*
            
            % sampleInputStream validity check
            if ~isa(sampleInputStream,'MultipleCompressedSampleInputStreamM')
                ME = MException('SpikeSaverM:IllegalConstructorArgumentException',...
                    ['Constructor argument sampleInputStream must be of java type ',...
                    'MultipleCompressedSampleInputStream.']);
                throw(ME);
            end
            
            % TODO Complete a path validity check for outputPath argument
            % savePath validity check
            if false
                ME = MException('SpikeSaverM:IllegalConstructorArgumentException',...
                    'Constructor argument outputPath is not a valid path.');
                throw(ME);
            else
                % Set the output path
                obj.savePath = savePath;
            end
            
            % Getting header from input stream
            header = sampleInputStream.getJavaHeader();
            
            % Create spikeFile object. outputPath should be corrected for
            % the argument should be a folder path and is used here as a file path
            
            % header does NOT always contain the correct packed array ID. See SpikeFinding.java
            
            % meanTimeConstant & spikeThreshold come from the argument HashMap
            % which likely comes from the xml files/user settings
            % These values should not be used on the way down by the saver except for storing
            % in header
            % So for now just set to 0
            meanTimeConstant = 0;
            spikeThreshold = 0;
            
            obj.spikeFile = SpikeFile([savePath, '.spikes'], header.getArrayID(),...
                meanTimeConstant, spikeThreshold,...
                header.getNumberOfSamples(),...
                header.getSamplingFrequency());
            
            obj.spikeSaver = SpikeSaver(obj.spikeFile);
            
        end
        
        function spikesFound = getSpikesFound(obj)
            spikesFound = obj.spikeSaver.getSpikesFound();
        end
        
        function processSpike(obj, spikeM)
            import edu.ucsc.neurobiology.vision.io.*
            
            if ~isa(spikeM,'SpikeM')
                ME = MException('SpikeSaverM:IllegalArgumentException',...
                    'processSpike argument spike is not a spikeM object.');
                throw(ME);
            end
            spike = Spike(spikeM.time,spikeM.electrode,spikeM.amplitude);
            
            obj.spikeSaver.processSpike(spike);
        end
        
        function finishSpikeProcessing(obj)
            % TODO deal with calling twice and avoid the java error.
           obj.spikeSaver.finishSpikeProcessing(); 
        end
    end
    
end

