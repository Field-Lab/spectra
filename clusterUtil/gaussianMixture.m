function [clusterIndexes, model, numClusters] = gaussianMixture( spikes )
    %GAUSSIANMIXTURE Clusters input data using variable k gaussian mixture model
    %
    %   Input:
    %       spikes: nSpikes x nDims array
    %
    %   Clustering with gaussian mixture, trying all values from 1 to maxGsn
    %   with maxGsn set in mVisionConfig.m
    %
    %   Finds best number of cluster with BIC thresholded metric
    %
    %   Output:
    %       clusterIndexes: nSpikes x 1 array, containing cluster index for
    %           each spike, with 0 for outliers if any.
    %       model: gmdistribution object, access methods provide cluster centroids,
    %           covariance matrices, ... Also allows reclustering of more data.
    %       numClusters: number of clusters for the data
    %
    % Author -- Vincent Deo -- Stanford University -- August 21, 2015
    
    %% Load mVision configuration
    config = mVisionConfig();
    clustConfig = config.getClustConfig();
    
    maxGsn = clustConfig.maxGaussians; % Maximum number of clusters
    
    %% Gaussian mixture, all k values clustering
    % Big data subsampling
    if size(spikes,1) > clustConfig.maxSpikes
        spikesToCluster = datasample(spikes,clustConfig.maxSpikes,1);
    else
        spikesToCluster = spikes;
    end
    
    GMmodels = cell(maxGsn,1); % Best clustering for varying number of clusters
    aic = zeros(1,maxGsn); % AIC metric
    bic = zeros(1,maxGsn); % BIC metric
    
    % Gaussian mixture fits
    for gsn = 1:maxGsn
        GMmodels{gsn} = fitgmdist(spikesToCluster,gsn,...
            'Options',statset('MaxIter',clustConfig.maxEMIter),...
            'Start','plus','RegularizationValue',clustConfig.regVal);
        aic(gsn) = GMmodels{gsn}.AIC;
        bic(gsn) = GMmodels{gsn}.BIC;
    end % gsn
    
    % Optimal number of cluster found
    % As first BIC value within 10% closer to final value reached for maximum
    % Number of clusters
    thr = 0.1;
    numClusters = find(bic < (1-thr)*bic(end)+thr*bic(1),1);
    
    model = GMmodels{numClusters}; % Final clustering model
    
    % Assigning clusters to spikes. clustConfig.clusterProb may be increased
    % to increase number of outliers and reduce clusters to their core
    clusterIndexes = (model.posterior(spikes) > clustConfig.clusterProb)*(1:numClusters)';
    
    %% Debugging plots
    if true
        %%
        % Plots 3D scatter plot of first 3 PCs, colored by cluster, and
        % along with kSig-sigma boxes of clusters. 
        figure(1)
        if size(spikes,1) > 10000
            dispsubset = randsample(size(spikes,1),10000);
        else
            dispsubset = 1:size(spikes,1);
        end
        scatter3(spikes(dispsubset,1),spikes(dispsubset,2),spikes(dispsubset,3),9,clusterIndexes(dispsubset));
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
        end % g
        numClusters;
    end % debug plots
    
end % function