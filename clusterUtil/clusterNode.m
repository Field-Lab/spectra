classdef clusterNode < handle
    %CLUSTERNODE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        low
        high
        numPoints
        
        pc @ pcPointer
        
        av
        covMat
        density
    end
    
    methods
        function obj = clusterNode(low, high, pc)
            validateattributes(low,{'numeric'},{'scalar'},'','lower bound',1);
            validateattributes(high,{'numeric'},{'scalar'},'','upper bound',2);
            validateattributes(high-low,{'numeric'},{'scalar','>=',0},'','cluster size');
            obj.low = low;
            obj.high = high;
            obj.numPoints = high - low + 1 ;
            
            obj.pc = pc;
            obj.update();
        end
        
        function boolVal = contains(obj, node)
            boolVal = obj.low <= node.low && obj.high >= node.high;
        end
        
        function val = overlap(obj,node)
            % This function check for an overlap between 2 clusters
            isOv = ~ ((obj.high < node.low) || (node.high < obj.low));
            % Measure overlap amount
            % Check case: inclusion
            if obj.contains(node) || node.contains(obj)
                val = true;
                return
            end
            % Non-inclusive overlap
            if isOv
                frac = 0.2; % Minimal fraction of overlap for either cluster to consider a merge 
                overLength = min(node.high - obj.low+1, obj.high-node.low+1);
                val = overLength / node.numPoints >= frac || overLength / obj.numPoints >= frac;
            else
                val = false;
            end
        end
        
        function newCluster = merge(obj,node)
            newCluster = clusterNode(min(obj.low,node.low),max(obj.high,node.high),obj.pc);
            newCluster.update();
        end
        
        function update(obj)
            obj.av =  mean(obj.pc.spikeComp(obj.pc.order(obj.low:obj.high),:),1);
            spkTemp = bsxfun(@minus,obj.pc.spikeComp(obj.pc.order(obj.low:obj.high),:),obj.av);
            obj.covMat =  spkTemp' * spkTemp / max(1,obj.numPoints-1);
            obj.density = (0.68^size(obj.pc.spikeComp,2))/prod(sqrt(eig(obj.covMat)));
        end
            
    end
    
end

