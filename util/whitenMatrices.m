function whitened = whitenMatrices( datapath, totSpikes, covMatrix, noiseCovMatrix )
    %WHITENMATRICES whitens the covariance matrix of spikes on an electrode
    % and neighbors using the associated covariance matrix of noise spikes
    % Based on vision's Noise-whitened-covariance
    %
    % Inputs:
    % covMatrix: nElectrodes x 1 cell array, each cell containing a
    %   nLPoints + nRPoints - 1 square matrix of covariance of aligned spikes
    % noiseCovMatrix: same, based on noise spike events.
    %
    % Output:
    % whitened: nElectrode x 1 cell array, each cell containing the whitened
    %   covariance matrix for the electrode
    %
    % Author -- Vincent Deo -- Stanford University -- August 11, 2015
    
    
    config = mVisionConfig(); covConfig = config.getCovConfig();
    
    % Not going to use data feed here, but need header properties
    dataSource = DataFileUpsampler(datapath);
    [adjacent,~] = catchAdjWJava( dataSource, covConfig.electrodeUsage);
    
    nElectrodes = dataSource.nElectrodes;
    
    % Adding zero event electrodes to disconnected list
    disconnected = or(totSpikes == 0, dataSource.disconnected);
    
    whitened = cell(nElectrodes,1);
    
    for el = 2:nElectrodes
        if disconnected(el)
            whitened{el} = zeros((covConfig.nRPoints+covConfig.nLPoints-1) * numel(adjacent{el}));
            continue
        end
        
        % Need to remove submatrices corresponding to zero-event electrodes, otherwise is not going
        % to work out fine.
        discardTags = repmat(disconnected(adjacent{el})',covConfig.nRPoints+covConfig.nLPoints-1,1);
        discardTags = discardTags(:)';
        
        % Whitening and replacement of disconnected blocks
        invMat = noiseCovMatrix{el}(~discardTags,~discardTags);
        if rank(invMat) < size(invMat,1)
            whitened{el}(~discardTags,~discardTags) = ...
                covMatrix{el}(~discardTags,~discardTags);
            continue
        end
        invMat = invMat^(-1/2);
        whitened{el} = zeros(size(covMatrix{el}));
        whitened{el}(~discardTags,~discardTags) = invMat * ...
            covMatrix{el}(~discardTags,~discardTags) * invMat;
    end % el
end

