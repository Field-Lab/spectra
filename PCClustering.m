function [clusterParams,neuronEls,neuronClusters,spikeTimesNeuron] = PCClustering(projSpikes, spikeTimesEl, varargin)
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
    
    nElectrodes = numel(projSpikes);
    clusterParams = cell(nElectrodes,1);
    
    
    neuronEls = cell(nElectrodes,1);
    neuronClusters = cell(nElectrodes,1);
    spikeTimesNeuron = cell(nElectrodes,1);
    
    %%%
    if nargin == 3
        %    load(varargin{1}) % reload model mode
    end
    %%%
    
%     tmp = neuronViewer('X:\EJgroup_data\NeurRawTest\vision\');
    
    % Optional parfor here - don't put if already parallelizing on files - put if single file processing
    % parfor
    for el = 2:nElectrodes
        if numel(projSpikes{el}) == 0
            continue
        end
        if size(projSpikes{el},2) < nDims
            throw(MException('',['PCClustering: insufficient dimensions (',num2str(size(projSpikes{el},2)),')saved in .prj.mat - please recompute projections.']));
        end
        
        try
            %%
            glmTimer = tic;
            
            %% Clustering function
            
            % Whitening
            [v,d] = eig(1/(size(projSpikes{el},1)) *...
                projSpikes{el}(:,1:nDims)' * projSpikes{el}(:,1:nDims));
            
            
            %             [clusterIndexes, model, numClusters] = gaussianMixture(projSpikes{el}(:,1:nDims) *...
            %                 v * diag(diag(d).^-0.5) * v');
            %             [clusterIndexes, model, numClusters] = spectralClustering(projSpikes{el}(:,1:nDims) *...
            %                 v * diag(diag(d).^-0.5) * v');
            
%             figure(4);
%             [colorRef,spTimesRef] = tmp.drawClustersColor(el-1,15000,0);
            
            [clusterIndexes, model, numClusters] = spectralClustering(projSpikes{el}(:,1:nDims));
            
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
                disp(sprintf(['Electrode ',num2str(el),':\n',...
                    num2str(numClusters),' neurons found over ',num2str(size(projSpikes{el},1)),' spikes\n',...
                    'Time for Electrode Clustering ',num2str(toc(glmTimer)),' seconds\n',...
                    '-----------------------']));
            end % if debug
            
        catch error
            disp(error);
            disp(['Error at electrode ',num2str(el),', skipping.']);
        end
        
    end % el
    
    % Concatenating all neurons found
    neuronEls = vertcat(neuronEls{:});
    neuronClusters = vertcat(neuronClusters{:});
    spikeTimesNeuron = vertcat(spikeTimesNeuron{:});
    
end % function
