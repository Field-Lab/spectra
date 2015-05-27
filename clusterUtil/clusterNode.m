classdef clusterNode < handle
    %CLUSTERNODE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        low
        high
        numPoints
        av
        sigma
        density
    end
    
    methods
        function obj = clusterNode(low, high, av, sigma)
            validateattributes(low,{'numeric'},{'scalar'},'','lower bound',1);
            validateattributes(high,{'numeric'},{'scalar'},'','upper bound',2);
            validateattributes(high-low,{'numeric'},{'scalar','>=',0},'','cluster size');
            obj.low = low;
            obj.high = high;
            obj.numPoints = high - low + 1 ;
            obj.av = av;
            obj.sigma = sigma;
            obj.density = (0.68^3)/prod(sigma);
        end
        
        function boolVal = contains(obj, node)
            boolVal = obj.low <= node.low && obj.high >= node.high;
        end
    end
    
end

