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
    
    %%
    spikes = double(spikes);
    %%
    
    % Determine if subsampling the data is necessary
    if specConfig.nSpikesL < size(spikes,1)
        subsampleTag = specConfig.lapAvIter;
    else
        subsampleTag = 1;
    end
    
    % sigma and max distance threshold of gaussian/inverse square metrics
    % are given as constant provided data is normalized so that the 1-sigma
    % box has n-dimensional volume 1.
    boxVolNorm = prod(std(spikes,1)) .^ (1/size(spikes,2)); % 1-sigma volume of data
    sigmaNorm = (specConfig.sigmaDist * boxVolNorm)^ 2; % Normalized sigma
    dthr = (specConfig.maxDistance * boxVolNorm ); % Normalized threshold
    
    % Define affinity metric 
    % Gaussian
    thr = exp(-dthr^2./(2*sigmaNorm));
    
    metric = @(euclDist2Norm) exp(-euclDist2Norm);
    % metric = @(euclDist2) euclDist2 <= dthr^2;
    % Inverse square
%     thr = sigmaNorm / (dthr.^2 + sigmaNorm);
%     metric = @(euclDist2) max(0, sigmaNorm./(euclDist2 + sigmaNorm) - thr) ./ (1-thr);
     
    % eigenvalues storage preallocation
    eigStack = zeros(subsampleTag,specConfig.subspaceDim);
    
    % Tags and counter for outlier criterion
    failTag = true;
    failCounter = 0;
    
    while failTag
        
        % Iterations to generate smoothed spectrum of graph laplacian
        % The eigenvalue gap is too unstable for a single data subset
        for iterIndex = 1:subsampleTag            
            if subsampleTag == 1 % data is small, not subsampling
                spikesToCluster = spikes;
            else % Subsampling
                if iterIndex < subsampleTag
                    spikesToCluster = datasample(spikes,specConfig.nSpikesL,1);
                else % At last iteration, pick larger subset for k-means precision
                    if specConfig.nSpikesK < size(spikes,1)
                        spikesToCluster = datasample(spikes,specConfig.nSpikesK,1);
                    else
                        spikesToCluster = spikes;
                    end
                end % iterIndex < subsampleTag
            end % subsampleTag == 1
            
            % Affinity matrix
            wmat = sum(bsxfun(@minus, permute(spikesToCluster,[1 3 2]),...
                permute(spikesToCluster,[3 1 2])).^2,3);
            
            sigmas = 1./sqrt(UtilSpectral.nthSmallest(0.03 * size(wmat,1),wmat));
%             sigmas = sqrt(UtilSpectral.nthSmallest(7,wmat));
            
            
%             sigprod = bsxfun(@times,sigmas',sigmas);
%             wmat = sigprod./(wmat+sigprod);
            wmat = exp(-bsxfun(@times,bsxfun(@times,wmat,sigmas'),sigmas));
            
            
            % Row weights
            dvector = sum(wmat,2);
            dinv = dvector .^ (-0.5);
            
            % Normalized graph (Id - Laplacian)
            % laplacian = eye(numel(dinv)) - bsxfun(@times,dinv',bsxfun(@times,dinv,wmat));
            laplacian = bsxfun(@times,dinv',bsxfun(@times,dinv,wmat));
            % laplacian = eye(numel(dinv)) - bsxfun(@times,dinv',wmat);
            
            % Extraction of dominant eigenvalues and eigenvectors
            options.issym = true; options.isreal = true; options.maxit = 1000;
            [evectors, evalues] = eigs(laplacian,specConfig.subspaceDim,'lm',options);
            evalues = 1 - evalues;
            
            [evalues,perm] = sort(diag(evalues)); % Sometimes, they come unsorted?
            evectors = evectors(:,perm);
            
            eigStack(iterIndex,:) = real(evalues); % Sometimes, eigs has imaginary residue?
        end % iterIndex
        
        % Smooth out laplacian spectrum
        sortedEigs = mean(eigStack,1);
        
        % Get spectrum gap index - spectral clustering best number of clusters
        gaps = sortedEigs(2:end)-sortedEigs(1:(end-1));
        [~,numClusters] = max(gaps);
        
        % Vision compatibily and avoid loophole of bad laplacian giving >> 10
        % number of clusters
        if numClusters > 15
            numClusters = 15;
        end
        
        % Dimensionality reduction and normalization of eigenvectors component by component
        evectorsNorm = evectors(:,1:min(size(evectors,2),numClusters));
        evectorsNorm = bsxfun(@rdivide,evectorsNorm,sqrt(sum(evectorsNorm.^2,2)));
        
        % Computing affinities of ALL points relative to subset used for laplacian
        % Then their PCs in Laplacian eigenspace
        % Proceeding by chunks of 1000 (RAM consuming intermediate step)
        
        numPerBuff = 1000;
        PCs = cell(ceil(size(spikes,1)/numPerBuff),1);
        
        for buff = 1:(size(PCs,1)-1)
            PCs{buff} = exp(-bsxfun(@times,sigmas'.^2,sum(bsxfun(@minus,...
                permute(spikes((numPerBuff*(buff-1)+1):(numPerBuff*buff),:),[1 3 2]),...
                permute(spikesToCluster,[3 1 2])).^2,3))) * evectorsNorm;
            
%             PCs{buff} = bsxfun(@times,sigmas'.^2,1./bsxfun(@plus,sigmas'.^2,(sum(bsxfun(@minus,...
%                 permute(spikes((numPerBuff*(buff-1)+1):(numPerBuff*buff),:),[1 3 2]),...
%                 permute(spikesToCluster,[3 1 2])).^2,3)))) * evectorsNorm;
        end
        PCs{end} = exp(-bsxfun(@times,sigmas'.^2,sum(bsxfun(@minus,...
                permute(spikes((numPerBuff*(size(PCs,1)-1)+1):end,:),[1 3 2]),...
                permute(spikesToCluster,[3 1 2])).^2,3))) * evectorsNorm;
        
%         PCs{end} = bsxfun(@times,sigmas'.^2,1./bsxfun(@plus,sigmas'.^2,(sum(bsxfun(@minus,...
%                 permute(spikes((numPerBuff*(size(PCs,1)-1)+1):end,:),[1 3 2]),...
%                 permute(spikesToCluster,[3 1 2])).^2,3)))) * evectorsNorm;
            
        PCs = vertcat(PCs{:});
        
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
    
%     % Whitening PCnorm
%     [v,d] = eig(1/(size(PCNorm,1)) *...
%         (PCNorm' * PCNorm));
%     PCNorm = PCNorm * v * diag(diag(d).^-0.5) * v';
    
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
    
    if false % Debug plots
        %%
        figure(1) % Laplacian spectrum
        plot(sortedEigs(1:specConfig.subspaceDim),'-o');

        figure(1);hold off;plot(eigStack');hold on; plot(eigStack(end,:),'r-*'); plot(sortedEigs,'b-*');
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
        scatter3(PCNorm(:,1),PCNorm(:,2),PCNorm(:,3),9,clusterIn);
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
        plot(x,1000*metric(x.^2),'r');
        hold off
        sigmaNorm
        numClusters
    end % Debug plots
end % function