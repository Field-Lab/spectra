classdef pcPointer < handle
    %SPIKEPOINTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        spikeComp
        order
    end
    
    methods
        function obj = pcPointer(spikeComp,order)
            obj.spikeComp = spikeComp;
            obj.order = order;
        end
    end
    
end

