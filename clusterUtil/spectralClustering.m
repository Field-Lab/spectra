function [clusterIndexes, model, numClusters] = spectralClustering( spikes )
    %GAUSSIANMIXTURE Clusters input data using spectrum-gap spectral clustering 
    %
    %   Input data should be nSpikes x nDims
    %
    %   Output is clustering indexes for the input spikes,
    %   And a spectral Clustering object
    
    config = mVisionConfig();
    specConfig = config.getSpectralConfig();
        
    % subsampling, reassigning should be done IN the function - same for GMM clustering
    % add for iter = 1:specConfig.repeatTimes
    if specConfig.nSpikes < size(spikes,1)
        spikesToCluster = datasample(spikes,specConfig.nSpikes,1);
    else
        spikesToCluster = spikes;
    end
    
    specConfig.sigmaDist
    
    sigmaNorm =  specConfig.sigmaDist^ 2 ;
    dthr = specConfig.maxDistance ;
    
    thr = exp(-dthr^2./(2*sigmaNorm));
    
    metric = @(euclDist2) max(0, exp(-euclDist2./(2*sigmaNorm))-thr) ./ (1-thr);
    
    wmat = metric(sum(bsxfun(@minus, permute(spikesToCluster,[1 3 2]),...
        permute(spikesToCluster,[3 1 2])).^2,3));
    
    dvector = sum(wmat,2);
    
%     [dmat,i] = sort(dvector,'descend');
%     dmat = diag(dmat);
%     wmat = wmat(i,i);
%     imagesc(wmat(i,i)),axis image,colorbar;
    
%   unnormalized laplacian
%   laplacian = dmat - wmat;
    
    dinv = dvector .^ (-0.5);
    
    laplacian = eye(numel(dinv)) - bsxfun(@times,dinv',bsxfun(@times,dinv,wmat));
%     laplacian = eye(numel(dinv)) - bsxfun(@times,dinv',wmat);
    
    [evectors, evalues] = eig(laplacian);
    
    % Those 2 lines may be unnecessary
    [sortedEigs,perm] = sort(abs(diag(evalues)));
    evectors = evectors(:,perm);
    
    gaps = sortedEigs(2:end)-sortedEigs(1:(end-1));
    [~,eigIndex] = max(gaps);
    
    numClusters = eigIndex
    
    
%     skip = find(sortedEigs > specConfig.eigThreshold, 1, 'first')-1;
    
    
    % PCs = wmat * evectors(:,skip:min(size(evectors,2),...
    %    (skip+specConfig.subspaceDim-1)));
    
    evectorsNorm = evectors(:,1:min(size(evectors,2),specConfig.subspaceDim));
    evectorsNorm = bsxfun(@rdivide,evectorsNorm,sqrt(sum(evectorsNorm.^2,2)));
    
    
    PCs = metric(sum(bsxfun(@minus, permute(spikes,[1 3 2]),...
        permute(spikesToCluster,[3 1 2])).^2,3)) * evectorsNorm;
    
    
%     figure(1)
%     plot(sortedEigs(1:30),'+');
    
    PCNorm = bsxfun(@rdivide,PCs,sqrt(sum(PCs.^2,2)));
%     PCNorm = PCs;
    
    [clusterIndexes,model] = kmeans(PCNorm,numClusters,...
        'replicates',specConfig.kmeansRep,...
        'MaxIter',specConfig.maxIter);
    
%     figure(2)
%     scatter3(spikes(:,1),spikes(:,2),spikes(:,3),9,clusterIndexes);
%     xlabel('x'),ylabel('y'),zlabel('z');
%     colorbar;
%     
%     figure(3);
%     scatter3(PCNorm(:,5),PCNorm(:,4),PCNorm(:,3),9,clusterIndexes);
%     xlabel('x'),ylabel('y'),zlabel('z');
    
end