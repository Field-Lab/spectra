function [clusterParams,neuronEls,neuronClusters,spikeTimesNeuron] = PCClustering(prjFilePath, spikeTimesEl)
    %PCCLUSTERING Clusters spikes over each electrode, outputs (raw) neurons
    %
    % Inputs:
    %   prjFilePath: Path to the analysis folder in which to find the .prj.mat file
    %   spikeTimesEl: nElectrodes x 1 cell array; spikeTimesEl{el} is a nSpike x 1
    %       array containing the spike times on electrode el
    %
    % Outputs:
    %   clusterParams: nElectrodes x 1 cell array; clusterParams{el} is a nNeurons x 1
    %       array of object/parameters describing the clusters
    %
    %   neuronEls: nNeurons x 1 array containing the electrode number for neurons
    %   neuronClusters: nNeurons x 1 array containing the cluster number describing neuron i
    %
    % --> clusterParams{neuronEls(i)}(neuronClusters(i)) describes a neuron
    %
    %   spikeTimesNeurons: nNeurons x 1 cell array; spikeTimeNeurons{neuron} is a
    %       nSpikesOfNeuron x 1 array containing the spike times of the neuron
    %
    % Clustering is done in a specified clustering subfunction that returns cluster numbers for a
    % given clustering space. This function to reindex all in a global dataset indexing
    % 
    % Author -- Vincent Deo -- Stanford University -- August 27, 2015

    % Load config
    config = mVisionConfig();
    clustConfig = config.getClustConfig();
    
    nDims = clustConfig.nDims;
    
    % Allocate
    nElectrodes = numel(spikeTimesEl);
    clusterParams = cell(nElectrodes,1);
    
    neuronEls = cell(nElectrodes,1);
    neuronClusters = cell(nElectrodes,1);
    spikeTimesNeuron = cell(nElectrodes,1);
    
    % Optional parfor here - don't put if already parallelizing on files (parallelCaller level)
    % parfor
    parfor el = 2:nElectrodes
        
        % Subfunction loadprj is necessary to maintain parfor body transparency.
        projSpikes = loadPrj(prjFilePath,el);
        
        % Not enough spikes - do not cluster
        if size(projSpikes,1) <= clustConfig.minSpikes
            continue
        end
        
        % Not enough dimensions stored - throw error
        if size(projSpikes,2) < nDims
            throw(MException('',['PCClustering: insufficient dimensions (',num2str(size(projSpikes{el},2)),')saved in .prj.mat - please recompute projections.']));
        end
        
        try % Catch errors on failure electrodes. Meant to not crash parfor loops when processing multiple files - quite unnecessary with current stabler version
            glmTimer = tic;
            
            % Clustering function
            [clusterIndexes, model, numClusters] = spectralClustering(projSpikes(:,1:nDims));
            
            % Assignments
            clusterParams{el} = model;
   
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
            
            projSpikes = []; % Gain some RAM by deallocating. Cannot use clear in parfor, cannot trust the GC.
        
        catch error
            disp(error);
            disp(['Error at electrode ',num2str(el),', skipping.']);
            throw(error)
        end
        
    end % el
    
    % Concatenating all neurons found
    neuronEls = vertcat(neuronEls{:});
    neuronClusters = vertcat(neuronClusters{:});
    spikeTimesNeuron = vertcat(spikeTimesNeuron{:});
    
end % function