function [clusterIndexes, model, numClusters] = spectralClustering( spikes )
    %GAUSSIANMIXTURE Clusters input data using spectrum-gap spectral clustering
    %
    %   Inputs
    %   spikes: nSpikes x nDims array containing data to cluster in rows
    %
    %   Outputs
    %   clusterIndexes: cluster number for each data point, in the range
    %       0 (outlier) - numClusters
    %   model: number of clusters x number of laplacian PCs matrix, containing
    %       cluster centroid coordinates in spectral space
    %   numclusters: resulting number of clusters
    %
    % Author -- Vincent Deo -- Stanford University -- August 5, 2015
    
    % Load mVision config and spetral custering configuration structure
    config = mVisionConfig();
    specConfig = config.getSpectralConfig();
    
    if ~strcmp(class(spikes),'double')
        spikes = double(spikes);
    end
    
    % sigma and max distance threshold of gaussian/inverse square metrics
    % are given as constant provided data is normalized so that the 1-sigma
    % box has n-dimensional volume 1.
%     boxVolNorm = prod(std(spikes,1)) .^ (1/size(spikes,2)); % 1-sigma volume of data
%     sigmaNorm = (specConfig.sigmaDist * boxVolNorm)^ 2; % Normalized sigma
%     dthr = (specConfig.maxDistance * boxVolNorm ); % Normalized threshold
    
    thrVal = 0.05;
    nNeighbors = 10;
    
    % Tags and counter for outlier criterion
    failTag = true;
    failCounter = 0;
    
    while failTag
        if specConfig.nSpikesK < size(spikes,1)
            spikesToCluster = datasample(spikes,specConfig.nSpikesK,1);
        else
            spikesToCluster = spikes;
        end
        
        % Distance matrix
        wmat = sum(bsxfun(@minus, permute(spikesToCluster,[1 3 2]),...
            permute(spikesToCluster,[3 1 2])).^2,3);
        
        
        %%%
        % Fixed or variable Kth neighbor ? multiplier factor ?
        sigmas = 1./sqrt(UtilSpectral.nthSmallest(nNeighbors,wmat)); % 0.03 * size(wmat,1)
%         sigmas = sqrt(UtilSpectral.nthSmallest(7,wmat));
        %%%
            
%         sigprod = bsxfun(@times,sigmas',sigmas);
%         wmat = sigprod./(wmat+sigprod);
        
        % Locally scaled affinity matrix
        % wmat = exp(-bsxfun(@times,bsxfun(@times,wmat,sigmas'),sigmas));
        wmat = 1./(1+bsxfun(@times,bsxfun(@times,wmat,sigmas'),sigmas).^2);
%         wmat(wmat < thrVal) = 0;
        
        % Row weights
        dvector = sum(wmat,2);
        dinv = dvector .^ (-0.5);
            
        % Normalized graph (Id - Laplacian)
        % laplacian = eye(numel(dinv)) - bsxfun(@times,dinv',bsxfun(@times,dinv,wmat));
        laplacian = bsxfun(@times,dinv',bsxfun(@times,dinv,wmat));
        % laplacian = eye(numel(dinv)) - bsxfun(@times,dinv',wmat);
        
        
        % Extraction of dominant eigenvalues and eigenvectors
        options.issym = true; options.isreal = true; options.maxit = 1000;
        [evectors, evalues] = eigs(laplacian,specConfig.maxClust + 5,'lm',options);
        evalues = 1 - evalues;
            
        [sortedEigs,perm] = sort(real(diag(evalues))); % Sometimes, they come unsorted, and with imaginary remainder
        evectors = evectors(:,perm);
        
        
        % Number of clusters determination
        currVectors = evectors(:,1:2); % Starting at 2 clusters
        quality = zeros(specConfig.maxClust-1,1);
        alignedVectors = cell(specConfig.maxClust-1,1);
        for C = 2:specConfig.maxClust
            [~,quality(C-1),alignedVectors{C-1}] = evrot(currVectors,1);
            currVectors = [alignedVectors{C-1},evectors(:,C+1)];
        end
        numClusters = find((max(quality)-quality) < 0.005,1,'last') + 1;
        
        % Dimensionality reduction and normalization of eigenvectors component by component
        % evectorsNorm = bsxfun(@rdivide,evectorsNorm,sqrt(sum(evectorsNorm.^2,2)));
        
        % Computing affinities of ALL points relative to subset used for laplacian
        % Then their PCs in Laplacian eigenspace
        % Proceeding by chunks of 1000 (RAM consuming intermediate step)
        
        evectors = evectors(:,1:numClusters);
        
        numPerBuff = 1000;
        PCs = cell(ceil(size(spikes,1)/numPerBuff),1);
        tic
        for buff = 1:(size(PCs,1)-1)
            adj = sum(bsxfun(@minus,...
                permute(spikes((numPerBuff*(buff-1)+1):(numPerBuff*buff),:),[1 3 2]),...
                permute(spikesToCluster,[3 1 2])).^2,3);
            sigmaBuff = 1./sqrt(UtilSpectral.nthSmallest(nNeighbors,adj)); % 0.03 * size(wmat,1)
            
            % PCs{buff} = exp(-bsxfun(@times,sigmaBuff,bsxfun(@times,sigmas',adj))) * evectors;
            adj = 1./(1+bsxfun(@times,sigmaBuff,bsxfun(@times,sigmas',adj)).^2);
%             adj(adj < thrVal) = 0;
            PCs{buff} = adj * evectors;
        end
        adj = sum(bsxfun(@minus,...
                permute(spikes((numPerBuff*(size(PCs,1)-1)+1):end,:),[1 3 2]),...
                permute(spikesToCluster,[3 1 2])).^2,3);
        
        sigmaBuff = 1./sqrt(UtilSpectral.nthSmallest(10,adj)); % 0.03 * size(wmat,1)
        adj = 1./(1+bsxfun(@times,sigmaBuff,bsxfun(@times,sigmas',adj)).^2);
%         adj(adj < thrVal) = 0;
        PCs{end} = adj * evectors;
        % PCs{end} = exp(-bsxfun(@times,sigmaBuff,bsxfun(@times,sigmas',adj))) * evectors;
        
        clear adj
        
        PCs = vertcat(PCs{:});
        toc
        
%         tic
%         PCs = UtilSpectral.computeSpectralProjections(nNeighbors,thrVal,spikesToCluster,sigmas,spikes,evectors);
%         toc
        
        % Normalize PCs on unit hypersphere
        % Discard outliers (unconnected to any laplacian-building point, hence
        % in kernel of eigenspace projector)
        PCNorm = bsxfun(@rdivide,PCs,sqrt(sum(PCs.^2,2)));
        discard = isnan(PCNorm(:,1));
        PCNorm(discard,:) = [];
        
        % Check for outlier fraction
        % If too high, then subset is unsatisfying, start over.
        if nnz(discard) > 0.01*size(PCs,1);
            if specConfig.debug
                disp('Failure at 1% max outlier criterion.');
            end
            failCounter = failCounter + 1;
        else
            failTag = false;
        end % nnz(discard) > 0.01*size(PCs,1)
        
        if failCounter >= 10
            disp('SpectralClustering: Could not meet 1% outlier rule.');
            failTag = false;
        end
    end % while failTag
    
    
    % K-means
    if size(PCNorm,1) > specConfig.maxPts
        % Dataset subsampling for k-means
        [PCNormRed,~] = datasample(PCNorm,specConfig.maxPts,1);
        
        [~,model] = kmeans(PCNormRed,numClusters,...
            'replicates',specConfig.kmeansRep,...
            'MaxIter',specConfig.maxIter);
        % Assigning nearest centroid to all points (even not used in k-means)
        [~,clusterIn] = min(sum(bsxfun(@minus, permute(PCNorm,[1 3 2]),...
                permute(model,[3 1 2])).^2,3),[],2);       
    else % Computing k-means directly on all points
        [clusterIn,model] = kmeans(PCNorm,numClusters,...
            'replicates',specConfig.kmeansRep,...
            'MaxIter',specConfig.maxIter);
    end
    
    % Assigning cluster 0 to outliers
    clusterIndexes = zeros(size(PCs,1),1);
    clusterIndexes(~discard) = clusterIn;
    
    if true % Debug plots
        %%
        figure(1) % Laplacian spectrum
        plot(quality,'+-');
        % plot(sortedEigs,'-o');

        %%
        figure(2) % Clusters in spike PC space
        if size(spikes,1) > 10000
            dispsubset = randsample(size(spikes,1),10000);
        else
            dispsubset = 1:size(spikes,1);
        end
        scatter3(spikes(dispsubset,1),spikes(dispsubset,2),spikes(dispsubset,3),9,clusterIndexes(dispsubset));
        xlabel('x'),ylabel('y'),zlabel('z');
        f = gca;
        f.DataAspectRatio = [1 1 1];
        caxis([-1 max(clusterIndexes)]); colormap hsv; colorbar;
        title('Spectral');
        %%
        figure(3); % Clusters in laplacian PC space.
        % Not very meaningful as PCs beyond 3 are still very significative
        if size(PCNorm,2) >= 3
            scatter3(PCNorm(:,1),PCNorm(:,2),PCNorm(:,3),9,clusterIn);
        else
            scatter(PCNorm(:,1),PCNorm(:,2),9,clusterIn);
        end
        xlabel('x'),ylabel('y'),zlabel('z');
        caxis([-1 max(clusterIndexes)]); colormap hsv; colorbar;
        title('Laplacian eigenspace');
        %%
        figure(4)
        tmat = sqrt(sum(bsxfun(@minus, permute(spikesToCluster,[1 3 2]),...
                permute(spikesToCluster,[3 1 2])).^2,3));
        hist(tmat(:),1000);
        hold on;
        x = min(tmat(:)):0.01:max(tmat(:));
%         plot(x,1000*metric(x.^2),'r');
        hold off
        numClusters
    end % Debug plots
end % function