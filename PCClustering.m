function [clusterParams,neuronEls,neuronClusters,spikeTimesNeuron] = PCClustering(projSpikes, spikeTimesEl)
%PCCLUSTERING Outputs dummy clusters
%
% Specifications
% Inputs
% projSpikes is a nElectrodes x 1 cell array
% projSpikes{el} is a nSpikes x nDims array containing the PC components of all spikes
%
% spikeTimesEl is a nElectrodes x 1 cell array
% spikeTimesEl{el} is a nSpike x 1 array containing the spike times on electrode el
%
% Outputs
% clusterParams is a nElectrodes x 1 cell array
% clusterParams{el} is either 1 object or a nClusters x 1 array of object/parameters describing the clusters
%
% neuronEls is a nNeurons x 1 array containing the electrode number for neurons
% neuronClusters is a nNeurons x 1 array containing the cluster number describing neuron i
% --> clusterParams{neuronEls(i)}(neuronClusters(i)) describes a neuron
%
% spikeTimesNeurons is a nNeurons x 1 cell array
% spikeTimeNeurons{neuron} is a nSpikesOfNeuron x 1 array containing the spike times of the
% neuron
%
% The clustering structure will allow to determine if any spike is the member of any cluster
% Then apply the requiring permutation/concatenation/discarding to build spikeTimeNeurons

dims = 3;

nElectrodes = numel(projSpikes);
clusterParams = cell(nElectrodes,1);
neuronEls = [];
neuronClusters = [];
spikeTimesNeuron = [];

for el = 2:nElectrodes
    try
        if numel(projSpikes{el}) == 0
            continue
        end
        %% Quick density computation by n-D binning
        dimMax = max(projSpikes{el},[],1);
        dimMin = min(projSpikes{el},[],1);
        dimBins = 30;
        binSpace = (dimMax-dimMin)/dimBins;
        binned = floor(bsxfun(@rdivide,projSpikes{el},binSpace));
        occupiedBins = size(unique(binned(:,1:dims),'rows'),1);
        
        avDens = size(projSpikes{el},1)./occupiedBins;
        
        
        %% OPTICS algorithm does the pre-clustering
        tic
        [ SetOfClusters, RD, CD, order ] = cluster_optics(projSpikes{el}(:,1:dims), round(5*avDens), 0);
        disp(['Time for optics pre-clustering ',num2str(toc)]);
        
        pc = pcPointer(projSpikes{el}(:,1:dims),order);
        
        store = [];
        % linear = zeros(1,numel(order));
        n = size(SetOfClusters,2);
        for k = 1:n
            store = [store,[SetOfClusters(k).start;SetOfClusters(k).end;SetOfClusters(k).end-SetOfClusters(k).start+1]];
        end
        
        %% OPTICS hierarchical clusters post-processing
        tree = clusterTree(store(1,:),store(2,:),pc);
        tree.reduce();
        
        [nodeArray,~] = tree.enumLeaves();
        
        relativeAvs = zeros(numel(nodeArray));
        
        for clusterIndex = 1:numel(nodeArray)
            cluster = nodeArray(clusterIndex);
            [v,d] = eig(cluster.covMat);
            for cl = 1:numel(nodeArray)
                relativeAvs(clusterIndex,cl) = norm(v' * ((nodeArray(cl).av'-cluster.av')./ sqrt(diag(d))));
            end
        end
        
        %% Discarding intersecting clusters
        relativeAvs(eye(numel(nodeArray)) == 1) = Inf;
        discard = false(1, numel(nodeArray));
        
        for cl = 1:numel(nodeArray)
            discard(cl) = any(relativeAvs(cl,:) <= 2);
            relativeAvs(cl:end,discard) = Inf;
        end
        
        nodeArrayRed = nodeArray(~discard);
        
        %% Gaussian mixture model
        gsn = numel(nodeArrayRed);
        if gsn > 8 % Vision fuckup to plan if we go around this limit, which is set in config.xml
            gsn = 8;
        end
        S.mu = zeros(gsn,dims);
        S.Sigma = zeros(dims,dims,gsn);
        S.ComponentProportion = zeros(1,gsn);
        
        for clusterIndex = 1:gsn
            cluster = nodeArrayRed(clusterIndex);
            S.mu(clusterIndex,:) = cluster.av';
            S.Sigma(:,:,clusterIndex) = cluster.covMat;
            S.ComponentProportion(clusterIndex) = cluster.numPoints./numel(order);
        end
        
        el
        %  R2015
        tic
        clusterParams{el} = fitgmdist(projSpikes{el}(:,1:dims),gsn,'Start',S,'RegularizationValue',0.01);
        %  R2014
        % clusterParams{el} = fitgmdist(projSpikes{el}(:,1:dims),gsn,'Start',S,'Regularize',0.01);
        disp(['Time for Gaussian Mixture Clustering ',num2str(toc)]);
        
        %% Assigning output
        neuronEls = [neuronEls;el*ones(gsn,1)];
        neuronClusters = [neuronClusters;(1:gsn)'];
        
        spikeClust = clusterParams{el}.posterior(projSpikes{el}(:,1:dims)) >= 0.8;
        temp = cell(gsn,1);
        for i = 1:gsn
            temp{i} = spikeTimesEl{el}(find(spikeClust(:,i)));
        end
        spikeTimesNeuron = [spikeTimesNeuron;temp];
    catch error
        disp(['Error at electrode ',el,', skipping.']);
    end
end % el
end % function

