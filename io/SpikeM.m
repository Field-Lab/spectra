classdef SpikeM < handle
    %SPIKEM Matlab class to encapsulate Spikes
    %   Built from scratch without encapsulating the java Spike class
    %   Follows the java class spike Template
    %   Matlab methods encapsulating java methods with argument a java Spike
    % should accept a SpikeM as argument then put it into a Spike object
    % before calling the java method
    
    properties (SetAccess = immutable, GetAccess = public)
        time        % time of occurence in samples
        electrode   % electrode number
        amplitude   % spike amplitude
    end
    
    methods
        function obj = SpikeM(time,electrode,amplitude)
            if nargin ~= 0
            validateattributes(time,{'numeric'},{'nonempty','integer','2d'},'','time',1);
            validateattributes(electrode,{'numeric'},{'nonempty','integer','size',size(time)},'','electrode',2);
            validateattributes(amplitude,{'numeric'},{'nonempty','size',size(time)},'','amplitude',3);
              [n,p] = size(time);
              obj(n,p) = SpikeM;
                for i = n:-1:1
                    for j = p:-1:1
                        obj(i,j).time = time(i,j);
                        obj(i,j).electrode = electrode(i,j);
                        obj(i,j).amplitude = amplitude(i,j);
                    end % j
                end % i
            end % nargin
        end % constructor
        
    end % methods
    
end

