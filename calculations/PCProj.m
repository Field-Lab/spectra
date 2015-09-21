function [projSpikes,eigenValues,eigenVectors,spikeTimes] = PCProj(dataPath, timeCommand, spikesTotal, covMatrix, averages, totSpikes)
    %PCPROJ digonalizes electrode covariance matrices, computes spikes principal components
    % 
    % Inputs:
    %   dataPath: path to data source folder
    %   timeCommand: time window of recordeddata to analyze 
    %   spikesTotal: nSpikes x 2 array; all dataset found spikes, as stored in .spikes.mat file.
    %       First column is spike time (sorted order) and second column is spike electrode
    %   covMatrix: nElectrodes x 1 cell array; covMatrix{el} contains the (whitened) covariance
    %       matrix of spikes on electrode el, of side nPoints*nNeighbors
    %   averages: nElectrodes x 1 cell array; averages{el} is a 1 x (nPoints*nNeighbors) array,
    %       containing average of all spikes on electrode el.
    %   totSpikes: nElectrodes x 1 array, containing total spike count per electrode
    %
    % Outputs:
    %   projSpikes: nElectrodes x 1 cell array; projSpikes{el} is a totSpikes(el) x nDims array
    %       containing spikes principal components
    %   eigenValues: nElectrodes x 1 cell array; eigenValues{el} is a nDims x 1 array containing the
    %       dominant eigenvalues for electrode el.
    %   eigenVectors: nElectrodes x 1 cell array; eigenValues{el} is a (nPoints*nNeighbors) x nDims array containing the
    %       dominant eigenvectors for electrode el.
    %   spikeTimes: nElectrodes x 1 cell array; spikeTimes{el} is a totSpikes(el) x 1 array containing (sorted) spike times on electrode el
    %       To be used with projSpikes{el}, row order matching.
    %
    % Author -- Vincent Deo -- Stanford University -- August 27, 2015
    
    
    %% Argument validation
    if ~(exist(dataPath,'file') == 2 || exist(dataPath,'file') == 7)
        throw(MException('','ProjCalculations: data folder|file does not exist'));
    end
    
    %% Load projections configuration
    config = mVisionConfig();
    projConfig = config.getProjConfig();
    
    rawDataSource = [dataPath,timeCommand];
    
    meanTimeConstant = projConfig.meanTimeConstant;
    
    nLPoints = projConfig.nLPoints;
    nRPoints = projConfig.nRPoints;
    
    nDims = projConfig.nDims;
    
    %% Creating data source
    dataSource = DataFileUpsampler(rawDataSource, meanTimeConstant, nLPoints, nRPoints);
    
    nElectrodes = dataSource.nElectrodes;
    disconnected = dataSource.disconnected;
    
    aligner = spikeAligner(dataSource); % Spike aligner for cubic spline resampling
    
    %% Projections storage and initialization
    eigenValues = cell(nElectrodes,1);
    eigenVectors = cell(nElectrodes,1);
    projSpikes = cell(nElectrodes,1);
    spikeTimes = cell(nElectrodes,1);
    currSpike = ones(nElectrodes,1);
    
    %% Electrode diagonalization loop
    for el = 2:nElectrodes
        if disconnected(el) % Skip case
            continue
        end
        
        % Diagonalize covariance matrix
        [v,d] = eig(covMatrix{el});
        [e,perm] = sort(diag(d),'descend');
        v = v(:,perm);
        eigenValues{el} = e(1:nDims);
        eigenVectors{el} = v(:,1:nDims);
        
        % Allocate cell arrays
        projSpikes{el} = zeros(totSpikes(el),nDims);
        spikeTimes{el} = zeros(1,totSpikes(el));
    end
    
    %% Main buffer loop
    while ~dataSource.isFinished % stopSample should be the first sample not loaded
        
        % Data source management
        [bufferStart,bufferEnd] = dataSource.loadNextBuffer();
        dataSource.createInterpolant();
        
        % Load Spikes
        spikesTemp = spikesTotal(and(spikesTotal(:,1) >= bufferStart + nLPoints, spikesTotal(:,1) < bufferEnd - nRPoints),:);
        
        % If no spikes at all are loaded, skip buffer
        if numel(spikesTemp) == 0
            continue
        end
        
        spikes = cell(nElectrodes,1);
        
        % Cell array-ify spikes
        for el = 1:nElectrodes
            spikes{el} = spikesTemp(spikesTemp(:,2) == el,1);
        end
        
        clear spikesTemp;
        
        
        %% Process by electrodes
        % Could parallel here, but actually slower due to the IO cost of sending to each worker.
        % The better part would be to change the while loop to a smarter parfor loop.
        for el = 2:nElectrodes
            
            nSpikes = numel(spikes{el});
            
            % Skip case
            if dataSource.disconnected(el) || nSpikes == 0;
                continue
            end

            % Store spike times
            spikeTimes{el}(currSpike(el):(currSpike(el) + nSpikes - 1)) = spikes{el};
            
            % Align all spikes with the aligner            
            spikesTemp = aligner.alignSpikes(el,spikes{el});
            
            % Compute PC projections
            projSpikes{el}(currSpike(el):(currSpike(el) + nSpikes - 1),:) = ...
                bsxfun(@minus,spikesTemp,averages{el}) * eigenVectors{el};
            
            % Increase spike counter
            currSpike(el) = currSpike(el) + nSpikes;
            
        end % el
    end % while ~isFinihed
    
end