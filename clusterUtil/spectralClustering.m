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
    if specConfig.nSpikesL < size(spikes,1)
        subsampleTag = specConfig.lapAvIter;
    else
        subsampleTag = 1;
    end
    
    specConfig.sigmaDist;
    
    boxVolNorm = prod(std(spikes,1)) .^ (1/size(spikes,2));
    sigmaNorm = (specConfig.sigmaDist * boxVolNorm)^ 2 ;
    dthr = (specConfig.maxDistance * boxVolNorm );
    
    thr = exp(-dthr^2./(2*sigmaNorm));
    
    metric = @(euclDist2) max(0, exp(-euclDist2./(2*sigmaNorm))-thr) ./ (1-thr);
    
    eigStack = zeros(subsampleTag,30);
    
    for iterIndex = 1:subsampleTag
        
        if subsampleTag == 1
            spikesToCluster = spikes;
        else
            if iterIndex < subsampleTag
                spikesToCluster = datasample(spikes,specConfig.nSpikesL,1);
            else
                if specConfig.nSpikesK < size(spikes,1)
                    spikesToCluster = datasample(spikes,specConfig.nSpikesK,1);
                else
                    spikesToCluster = spikes;
                end
            end
        end
        
        wmat = metric(sum(bsxfun(@minus, permute(spikesToCluster,[1 3 2]),...
            permute(spikesToCluster,[3 1 2])).^2,3));
        
        dvector = sum(wmat,2);
        %
        %     [dmat,i] = sort(dvector,'descend');
        %     dmat = diag(dmat);
        %     wmat = wmat(i,i);
        %     imagesc(wmat),axis image,colorbar;
        
        %   unnormalized laplacian
        %   laplacian = dmat - wmat;
        
        dinv = dvector .^ (-0.5);
        
%         laplacian = eye(numel(dinv)) - bsxfun(@times,dinv',bsxfun(@times,dinv,wmat));
        laplacian = bsxfun(@times,dinv',bsxfun(@times,dinv,wmat));
        %     laplacian = eye(numel(dinv)) - bsxfun(@times,dinv',wmat);
        
        %     [evectors, evalues] = eig(laplacian);
        options.issym = true;
        [evectors, evalues] = eigs(laplacian,30,'lm',options);
        evalues = 1 - evalues;
        
        [evalues,perm] = sort(diag(evalues));
        evectors = evectors(:,perm);
        % Those 2 lines may be unnecessary
        
        eigStack(iterIndex,:) = real(evalues);
        % iterIndex
    end
    
    sortedEigs = mean(eigStack,1);
    
    gaps = sortedEigs(2:end)-sortedEigs(1:(end-1));
    [~,eigIndex] = max(gaps);
    
    numClusters = eigIndex;
    %     skip = find(sortedEigs > specConfig.eigThreshold, 1, 'first')-1;
    
    
    % PCs = wmat * evectors(:,skip:min(size(evectors,2),...
    %    (skip+specConfig.subspaceDim-1)));
    
    evectorsNorm = evectors(:,1:min(size(evectors,2),specConfig.subspaceDim));
    evectorsNorm = bsxfun(@rdivide,evectorsNorm,sqrt(sum(evectorsNorm.^2,2)));
    
    
    PCs = metric(sum(bsxfun(@minus, permute(spikes,[1 3 2]),...
        permute(spikesToCluster,[3 1 2])).^2,3)) * evectorsNorm;
    
    if numClusters > 15
        numClusters = 15;
    end
    
    
    PCNorm = bsxfun(@rdivide,PCs,sqrt(sum(PCs.^2,2)));
    PCNorm(isnan(PCNorm(:,1)),:) = [];
    %     PCNorm = PCs;
    
    [clusterIndexes,model] = kmeans(PCNorm,numClusters,...
        'replicates',specConfig.kmeansRep,...
        'MaxIter',specConfig.maxIter);
    
    if false
        %%
        figure(1)
        plot(sortedEigs(1:30),'-o');
        
        figure(2)
        scatter3(spikes(:,1),spikes(:,2),spikes(:,3),9,clusterIndexes);
        xlabel('x'),ylabel('y'),zlabel('z');
        colorbar;
        
        figure(3);
        scatter3(PCNorm(:,1),PCNorm(:,2),PCNorm(:,3),9,clusterIndexes);
        xlabel('x'),ylabel('y'),zlabel('z');
    end
    
end