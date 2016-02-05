function [clusterIndexes, model, numClusters] = spectralClustering( spikes )
    %SPECTRALCLUSTERING Clusters input data using full-automated spectral clustering
    %
    %   Input:
    %       spikes: nSpikes x nDims array
    %
    %   Clustering with spectral clustering, with k-th neighbor distance setting a local
    %   sigma for each datapoint. Affinity matrix uses gaussian metrix over
    %   normalized euclidian distance.
    %   Finds best number of clusters with laplacian subspace eigenvector
    %   optimal rotations
    %
    %   Output:
    %       clusterIndexes: nSpikes x 1 array, containing cluster index for
    %           each spike, with 0 for outliers if any.
    %       model: k-means centroid matrix in spectral space. Not really useful
    %           out of this spectral framework.
    %       numClusters: number of clusters for the data
    %
    % Author -- Vincent Deo -- Stanford University -- August 21, 2015
    
    %% Load mVision configuration - prep data
    config = mVisionConfig();
    specConfig = config.getSpectralConfig();
    
    % Cast in double if not - eigs does not support single (R2015a)
    if ~isa(spikes,'double')
        spikes = double(spikes);
    end
    
    % Load java nth-neighbor distance computation loop
    % If PCClustering runs in parfor mode, there is some weird stuff happening
    % in sharing the java path with workers
    if ~exist('UtilSpectral','class')
        javaaddpath(['.',filesep,'clusterUtil']);
        if ~exist('UtilSpectral','class')
            throw(MException('','spectralClustering:java compatibility issues, check compile and runtime jvms.'));
        end
    end
    
    %% Big data subsampling for laplacian computation (quadratic in time and memory)
    if specConfig.spikesForLaplacian < size(spikes,1)
        spikesToCluster = datasample(spikes,specConfig.spikesForLaplacian,1);
    else
        spikesToCluster = spikes;
    end
    
    %% Computing and diagonalizing the graph Laplacian
    
    % Euclidian distance matrix
    wmat = sum(bsxfun(@minus, permute(spikesToCluster,[1 3 2]),...
        permute(spikesToCluster,[3 1 2])).^2,3);
    
    % Computing k-th neighbor sigma distance in java auxiliary loop
    sigmas = 1./sqrt(UtilSpectral.nthSmallest(specConfig.nthNeighbor,wmat));
    
    % Locally scaled affinity matrix
    wmat = exp(-bsxfun(@times,bsxfun(@times,wmat,sigmas'),sigmas));
    
    % Row weights
    dvector = sum(wmat,2);
    dinv = dvector .^ (-0.5);
    
    % Normalized (Id - Laplacian) of the graph
    laplacian = bsxfun(@times,dinv',bsxfun(@times,dinv,wmat)); % D^(-1/2) * A * D^(-1/2)
    laplacian = (laplacian + laplacian.')/2;
    
    % Extraction of dominant eigenvalues and eigenvectors
    options.issym = 1; options.isreal = 1;
    options.maxit = specConfig.eigsMaxIter;
    options.tol = specConfig.eigsTol;
    [evectors, evalues] = eigs(laplacian,specConfig.maxEV + 5,'lm',options);
    evalues = 1 - evalues; % Eigenvalues of Laplacian
    
    % Sorting and suppressing imaginary residue
    [~,perm] = sort(real(diag(evalues)));
    evectors = evectors(:,perm);
    
    %% Number of clusters determination
    % Done by iteratively finding optimal rotation of eigenvectors in
    % smallest eigenvalues spaces
    
    currVectors = evectors(:,1:2); % Starting at 2 clusters
    quality = zeros(specConfig.maxEV-1,1); % Array of quality metric
    alignedVectors = cell(specConfig.maxEV-1,1); % Rotated vectors at each iteration
    for C = 2:specConfig.maxEV
        [~,quality(C-1),alignedVectors{C-1}] = evrot(currVectors,1);
        currVectors = [alignedVectors{C-1},evectors(:,C+1)];
    end % C
    
    % Qualities are a 0-1 metrix, often max around 1
    % Best number of clusters is highest number within 0.02 of maximum
    numClusters = find((max(quality)-quality) < specConfig.qualityTol,1,'last') + 1;
    
    %% Dimensionality reduction - Laplacian PC space of dimension nClusters
    
    % evectors = evectors(:,1:numClusters); % Projection vectors
    evectors = alignedVectors{numClusters-1};
    
    % Projecting ALL spikes in the Laplacian eigenspace
    % Output array is nSpikes * nClusters
    % But intermediate affinity is of size nSpikes * size(laplacian,1)
    % Done in successive buffers to not blow up RAM
    
    numPerBuff = specConfig.reprojBufferSize; % spike buffer size
    PCs = cell(ceil(size(spikes,1)/numPerBuff),1); % PCs storage
    
    % Buffer loop
    for buff = 1:size(PCs,1)
        buffStart = numPerBuff*(buff-1)+1;
        buffEnd = min(size(spikes,1),numPerBuff*buff);
        
        adj = sum(bsxfun(@minus,...
            permute(spikes(buffStart:buffEnd,:),[1 3 2]),...
            permute(spikesToCluster,[3 1 2])).^2,3); % Adjacency matrix
        
        sigmaBuff = 1./sqrt(UtilSpectral.nthSmallest(specConfig.nthNeighbor,adj)); % Local scales
        
        adj = exp(-bsxfun(@times,sigmaBuff,bsxfun(@times,sigmas',adj))); % Affinity matrix
        % Add possible thresholding here
        PCs{buff} = bsxfun(@times,dinv',bsxfun(@times,sum(adj,2).^(-1/2),adj)) * evectors; % Projections
    end % buff
    
    clear adj
    
    PCs = vertcat(PCs{:}); % Compacting whole nSpikes * nClusters in one array
    
    %% Normalize PCs on unit hypersphere - discard outliers
    % Discard outliers (unconnected to any laplacian-building point,
    % hence in kernel of eigenspace projector)
    PCNorm = bsxfun(@rdivide,PCs,sqrt(sum(PCs.^2,2)));
    discard = isnan(PCNorm(:,1));
    PCNorm(discard,:) = [];
    
    % Check for outlier fraction
    % If too high, then subset is unsatisfying, start over.
    % note: cannot fail if no thresholding
    % Exception is handled in PCClustering.m,
    % which retries if either eigs or the outlier criterion crash.
    if nnz(discard) > specConfig.maxOutlierFraction*size(PCs,1);
        throw(MException('',...
            sprintf('Failure at %2.2f%% max outlier criterion.\n',...
            100*specConfig.maxOutlierFraction)));
    end
    
    %% k-means clustering in spectral space
    
    if size(PCNorm,1) > specConfig.spikesForKMeans % Too many spikes to cluster - subsample
        [PCNormRed,~] = datasample(PCNorm,specConfig.spikesForKMeans,1); % subset
        
        [~,modelSpectralSpace] = kmeans(PCNormRed,numClusters,...
            'replicates',specConfig.kMeansReplicas,...
            'MaxIter',specConfig.kMeansIterations); % k-means
        
        % Assigning nearest centroid to all points (even not used in k-means)
        [~,clusterIn] = min(sum(bsxfun(@minus, permute(PCNorm,[1 3 2]),...
            permute(modelSpectralSpace,[3 1 2])).^2,3),[],2);
    else % Computing k-means directly on all points
        [clusterIn,~] = kmeans(PCNorm,numClusters,...
            'replicates',specConfig.kMeansReplicas,...
            'MaxIter',specConfig.kMeansIterations);
    end
    
    % Assigning cluster 0 to outliers
    clusterIndexes = zeros(size(PCs,1),1);
    clusterIndexes(~discard) = clusterIn;
    
    model.nDimensions = size(spikes,2);
    model.numClusters = numClusters;
    model.centroids = zeros(numClusters,size(spikes,2));
    model.covariances = zeros(numClusters,size(spikes,2));
    model.mixFrac = zeros(numClusters,1);
    for c = 1:numClusters
        model.centroids(c,:) = mean(spikes(clusterIndexes == c,:),1);
        model.covariances(c,:) = var(spikes(clusterIndexes == c,:),1);
        model.mixFrac(c) = sum(clusterIndexes == c)./size(spikes,1);
    end
    
    %% Debug plots
    if false
        %% Quality function
        figure(1)
        plot(quality,'+-');
        
        %% Clusters in spikes PC space
        figure(2)
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
        %% Clusters in laplacian PC space.
        % Hard to read as what you see most of the time is a nClusters-dimensional
        % hypersphere projected in dimension 3
        figure(3);
        if size(PCNorm,2) >= 3
            scatter3(PCNorm(:,1),PCNorm(:,2),PCNorm(:,3),9,clusterIn);
        else
            scatter(PCNorm(:,1),PCNorm(:,2),9,clusterIn);
        end
        xlabel('x'),ylabel('y'),zlabel('z');
        caxis([-1 max(clusterIndexes)]); colormap hsv; colorbar;
        title('Laplacian eigenspace');
        %% Histogram of distance matrix
        figure(4)
        tmat = sqrt(sum(bsxfun(@minus, permute(spikesToCluster,[1 3 2]),...
            permute(spikesToCluster,[3 1 2])).^2,3));
        hist(tmat(:),1000);
        hold on;
        x = min(tmat(:)):0.01:max(tmat(:));
        hold off
        fprintf('%u clusters found\n',numClusters);
        fprintf('%u spikes processed\n',size(spikes,1));
    end % Debug plots
end % function
