classdef spikeAligner < handle
    %SPIKEALIGNER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        spikeBuffer
        
        dataSource
        
        nLPoints
        nRPoints
        nPoints
        
        adjacent
        maxAdjacent
        
        upSampRatio
        upSampStep
        
        resampleBase
    end
    
    methods
        function obj = spikeAligner(dataSource)
            validateattributes(dataSource,{'DataFileUpsampler'},{},'','dataSource',1);
            
            obj.dataSource = dataSource;
            
            config = mVisionConfig();
            alignConfig = config.getCovConfig();
            
            obj.nLPoints = alignConfig.nLPoints;
            obj.nRPoints = alignConfig.nRPoints;
            obj.nPoints = obj.nLPoints + obj.nRPoints + 1;
            
            [obj.adjacent,obj.maxAdjacent] = catchAdjWJava( dataSource, alignConfig.electrodeUsage);
            
            obj.upSampRatio = dataSource.upSampleRatio;
            obj.upSampStep = 1/obj.upSampRatio;
            
            obj.resampleBase = 0:obj.upSampStep:2;
            
            obj.spikeBuffer = zeros(1,(obj.nPoints-2) * obj.maxAdjacent);
            
        end % Constructor
        
        function alignedSpikes = alignSpikes(obj, el, spikes)
            
            nSpikes = numel(spikes);
            
            %% Check if buffer expansion is required
            if nSpikes > size(obj.spikeBuffer,1);
                obj.spikeBuffer = zeros(nSpikes,(obj.nPoints-2) * obj.maxAdjacent);
            end
            
            interpIndex = bsxfun(@plus,obj.resampleBase,double(spikes) -  obj.dataSource.bufferStart);
            interpSpikes = obj.dataSource.interpolant{el}(interpIndex(:));
            interpSpikes = reshape(interpSpikes,size(interpIndex));
            
            [~,offset] = min(interpSpikes,[],2);
            
            interpPoints = bsxfun(@plus,...
                (offset-1)/obj.upSampRatio + double(spikes) - obj.dataSource.bufferStart - obj.nLPoints,...
                1:(obj.nPoints-2));
            interpPointsLin = interpPoints(:);
            
            s = size(interpPoints);
            for elAdjIndex = 1:numel(obj.adjacent{el})
                elAdj = obj.adjacent{el}(elAdjIndex);
                obj.spikeBuffer(1:nSpikes,((obj.nPoints-2)*(elAdjIndex-1)+1):((obj.nPoints-2)*elAdjIndex)) =...
                    reshape(obj.dataSource.interpolant{elAdj}(interpPointsLin),s);
            end
            
            alignedSpikes = obj.spikeBuffer(1:nSpikes,1:(numel(obj.adjacent{el})*(obj.nPoints-2)));
            
            if false % Alignment debug plots
            %%
                clf
                for elAdjIndex = 1:numel(obj.adjacent{el})
                    elAdj = obj.adjacent{el}(elAdjIndex);
                    hold on
                    plot(1:size(obj.dataSource.rawData,2),obj.dataSource.rawData(elAdj,:)+(elAdjIndex-1)*150,'k+');
                    plot(1:0.05:size(obj.dataSource.rawData,2),...
                        obj.dataSource.interpolant{elAdj}(1:0.05:size(obj.dataSource.rawData,2))+(elAdjIndex-1)*150,'b--')
                    for sp = 1:numel(spikes)
                        plot(interpPoints(sp,:),...
                        alignedSpikes(sp,((elAdjIndex-1)*(obj.nPoints-2) + 1):(elAdjIndex*(obj.nPoints-2)))+(elAdjIndex-1)*150,'r+-');
                    end
                    offset
                    hold off
                end
            end
        end
    end
    
end

