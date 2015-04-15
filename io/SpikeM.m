classdef SpikeM
    %SPIKEM Matlab class to encapsulate Spikes
    %   Built from scratch without encapsulating the java Spike class
    %   Follows the java class spike Template
    %   Matlab methods encapsulating java methods with argument a java Spike
    % should accept a SpikeM as argument then put it into a Spike object
    % before calling the java method
    
    properties (SetAccess = immutable, GetAccess = public)
        time@int32        % time of occurence in samples
        electrode@int32   % electrode number
        amplitude@double  % spike amplitude
    end
    
    methods
        function obj = SpikeM(time,electrode,amplitude)
            obj.time = int32(time);
            obj.electrode = int32(electrode);
            obj.amplitude = double(amplitude);
        end        
    end
    
end

