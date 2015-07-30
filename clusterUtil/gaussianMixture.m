function [clusterIndexes, model, numClusters] = gaussianMixture( spikes )
    %GAUSSIANMIXTURE Clusters input data using variable k gaussian mixture model 
    %
    %   Input data should be nSpikes x nDims
    %   Clustering with gaussian mixture, trying all values from 1 to maxGsn
    %   And finding best one with a certain metric
    %
    %   Output is clustering indexes for the input spikes,
    %   And a gmdistribution object
    
    config = mVisionConfig();
    clustConfig = config.getClustConfig();
    
    maxGsn = clustConfig.maxGaussians;
    
    
    if size(spikes,1) > clustConfig.maxSpikes
        spikesToCluster = datasample(spikes,clustConfig.maxSpikes,1);
    else
        spikesToCluster = spikes;
    end
    
    
    GMmodels = cell(maxGsn,1);
    aic = zeros(1,maxGsn);
    bic = zeros(1,maxGsn);
    
    for gsn = 1:maxGsn
        GMmodels{gsn} = fitgmdist(spikesToCluster,gsn,...
            'Options',statset('MaxIter',clustConfig.maxEMIter),...
            'Start','plus','RegularizationValue',clustConfig.regVal);
        aic(gsn) = GMmodels{gsn}.AIC;
        bic(gsn) = GMmodels{gsn}.BIC;
    end
            
    thr = 0.1;
    numClusters = find(bic < (1-thr)*bic(end)+thr*bic(1),1);
           
    model = GMmodels{numClusters};
            
    %% Assigning output
            
    clusterIndexes = (model.posterior(spikes) > clustConfig.clusterProb)*(1:numClusters)';
    
end