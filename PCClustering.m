function [clusterParams,neuronEls,neuronClusters,spikeTimesNeuron] = PCClustering(prjFilePath, spikeTimesEl)
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
    
    config = mVisionConfig();
    clustConfig = config.getClustConfig();
    
    nDims = clustConfig.nDims;
    
    nElectrodes = numel(spikeTimesEl);
    clusterParams = cell(nElectrodes,1);
    
    
    neuronEls = cell(nElectrodes,1);
    neuronClusters = cell(nElectrodes,1);
    spikeTimesNeuron = cell(nElectrodes,1);
    
%     tmp = neuronViewer('X:\EJgroup_data\NeurRawTest\vision\');
    
    % Optional parfor here - don't put if already parallelizing on files - put if single file processing
    % parfor
    parfor el = 2:nElectrodes
        
        % Subfunction loadprj is necessary to maintain parfor body transparency.
        projSpikes = loadPrj(prjFilePath,el);
            
        if numel(projSpikes) <= clustConfig.minSpikes
            continue
        end
        if size(projSpikes,2) < nDims
            throw(MException('',['PCClustering: insufficient dimensions (',num2str(size(projSpikes{el},2)),')saved in .prj.mat - please recompute projections.']));
        end
        
%         try
            %%
            glmTimer = tic;
            
            %% Clustering function
            
            % Whitening
            %             [v,d] = eig(1/(size(projSpikes{el},1)) *...
            %                 projSpikes{el}(:,1:nDims)' * projSpikes{el}(:,1:nDims));
            %             [clusterIndexes, model, numClusters] = gaussianMixture(projSpikes{el}(:,1:nDims) *...
            %                 v * diag(diag(d).^-0.5) * v');
            %             [clusterIndexes, model, numClusters] = spectralClustering(projSpikes{el}(:,1:nDims) *...
            %                 v * diag(diag(d).^-0.5) * v');
            
            [clusterIndexes, model, numClusters] = spectralClustering(projSpikes(:,1:nDims));
            
            clusterParams{el} = model;
            
            
            % model should have a method assign - and we must do posterior assignments for spikes
            % not used in the clustering.
            
            neuronEls{el} = el*ones(numClusters,1);
            neuronClusters{el} = (1:numClusters)';
            
            spikeTimesNeuron{el} = cell(numClusters,1);
            
            for gsn = 1:numClusters
                spikeTimesNeuron{el}{gsn} = spikeTimesEl{el}(clusterIndexes == gsn);
            end
            
            %% Debug - plots
            if clustConfig.debug
                fprintf(['Electrode %u:\n%u neurons found over %u spikes\n',...
                    'Time for Electrode Clustering %f seconds\n',...
                    '-----------------------\n'],el,numClusters,size(projSpikes,1),toc(glmTimer));
            end % if debug
            
            projSpikes = [];
%         catch error
%             disp(error);
%             disp(['Error at electrode ',num2str(el),', skipping.']);
%             throw(error)
%         end
        
    end % el
    
    % Concatenating all neurons found
    neuronEls = vertcat(neuronEls{:});
    neuronClusters = vertcat(neuronClusters{:});
    spikeTimesNeuron = vertcat(spikeTimesNeuron{:});
    
end % function
