classdef spikeAligner < handle
    %SPIKEALIGNER Utility class - handles cubic spline subsample realignment of spikes
    % during covariance calculation and projections computation - provides standard behavior
    % 
    % Author -- Vincent Deo -- Stanford University -- August 27, 2015
    
    properties
        spikeBuffer % Allocated buffer
        
        dataSource % DataFileUpsampler, data feed
        
        % Alignment setup
        nLPoints
        nRPoints
        nPoints
        
        % electrodes setup
        adjacent
        maxAdjacent
        
        % Alignment precision
        upSampRatio
        upSampStep
        resampleBase
    end
    
    methods
        % Constructor
        %
        % Input:
        %   dataSource: DataFileUpsampler object, sample stream on which spikes will be aligned
        function obj = spikeAligner(dataSource)
            validateattributes(dataSource,{'DataFileUpsampler'},{},'','dataSource',1);
            
            obj.dataSource = dataSource;
            
            % Load config
            global GLOBAL_CONFIG
            alignConfig = GLOBAL_CONFIG.getCovConfig();
            
            % Set aligner config
            obj.nLPoints = alignConfig.nLPoints;
            obj.nRPoints = alignConfig.nRPoints;
            obj.nPoints = obj.nLPoints + obj.nRPoints + 1;
            
            % Load adjacency pattern
            [obj.adjacent,obj.maxAdjacent] = catchAdjWJava( dataSource, alignConfig.electrodeUsage);
            
            % Upsampling precision
            obj.upSampRatio = dataSource.upSampleRatio;
            obj.upSampStep = 1/obj.upSampRatio;
            
            obj.resampleBase = 0:obj.upSampStep:2;
            
            % Allocated buffer
            obj.spikeBuffer = zeros(1,(obj.nPoints-2) * obj.maxAdjacent);
            
        end % Constructor
        
        
        % Aligns the provided spikes on the available data with sub-sample precision
        % Performs cubic-spline interpolation and places the spike minimum exactly on a sample
        %
        % Input:
        %   el: electrode number - scalar
        %   spikes: list of spikes times on electrode el - column
        %
        % Return: 
        %   Aligned spikes: nSpikes x (nPoints*nNeighbor) array containing realigned spikes in rows
        function alignedSpikes = alignSpikes(obj, el, spikes)
            
            nSpikes = numel(spikes);
            
            %% Check if buffer expansion is required
            if nSpikes > size(obj.spikeBuffer,1);
                obj.spikeBuffer = zeros(nSpikes,(obj.nPoints-2) * obj.maxAdjacent);
            end
            
            % Compute interpolation indexes near the minimum, request cubic spline interpolation values at these
            % points to the data source
            interpIndex = bsxfun(@plus,obj.resampleBase,double(spikes) -  obj.dataSource.bufferStart);
            interpSpikes = obj.dataSource.interpolant{el}(interpIndex(:));
            interpSpikes = reshape(interpSpikes,size(interpIndex));
            
            % Find the minimum
            [~,offset] = min(interpSpikes,[],2);
            
            % Find 1-sample spaced realigned interpolation points along the spike
            interpPoints = bsxfun(@plus,...
                (offset-1)/obj.upSampRatio + double(spikes) - obj.dataSource.bufferStart - obj.nLPoints,...
                0:(obj.nPoints-3));
            interpPointsLin = interpPoints(:);
            
            s = size(interpPoints);
            
            % Successively load aligned electrode spikes and neighbor spikes in the buffer
            % All spikes performed at once, and put in successive, left to right (nSpikes x nPoints)
            % blocks of the buffer
            for elAdjIndex = 1:numel(obj.adjacent{el})
                elAdj = obj.adjacent{el}(elAdjIndex);
                obj.spikeBuffer(1:nSpikes,((obj.nPoints-2)*(elAdjIndex-1)+1):((obj.nPoints-2)*elAdjIndex)) =...
                    reshape(obj.dataSource.interpolant{elAdj}(interpPointsLin),s);
            end % elAdjIndex
            
            % Extract relevant lines of buffer
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
                end % elAdjIndex
            end % debug plots
        end % alignedSpikes
    end % methods
    
end % spikeAligner

