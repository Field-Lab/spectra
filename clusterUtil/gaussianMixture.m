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
    
    if true % Debug plots
        %%
        figure(4)
        scatter3(spikes(:,1),spikes(:,2),spikes(:,3),9,clusterIndexes);
        colormap jet
        colorbar
        title('Gaussian mixture')
        
        for g = 1:numClusters
            covMat = model.Sigma(1:3,1:3,g);
            center = model.mu(g,1:3);
            [v,d] = eig(covMat);
            
            kSig = 1;
            lineMat = [1,-1,-1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,-1,1;...
                -1,-1,1,1,1,1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1;...
                -1,-1,-1,-1,1,1,1,1,-1,-1,1,-1,1,1,-1,1,-1];
            cube = bsxfun(@plus,center',v * kSig * bsxfun(@times,sqrt(diag(d)),lineMat));
            
            hold on
            plot3(cube(1,:),cube(2,:),cube(3,:),'k-','linewidth',2);
            axis tight
            hold off
        end
        numClusters;
    end % debug plots
    
end