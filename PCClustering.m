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
    
    config = mVisionConfig();
    clustConfig = config.getClustConfig();
    
    nDims = clustConfig.nDims;
    maxGsn = clustConfig.maxGaussians;
    
    nElectrodes = numel(projSpikes);
    clusterParams = cell(nElectrodes,1);
    
    
    neuronEls = cell(nElectrodes,1);
    neuronClusters = cell(nElectrodes,1);
    spikeTimesNeuron = cell(nElectrodes,1);
    
    parfor el = 2:nElectrodes
        if numel(projSpikes{el}) == 0
            continue
        end
        try
            %%
            glmTimer = tic;
            
            GMmodels = cell(maxGsn,1);
            aic = zeros(1,maxGsn);
            bic = zeros(1,maxGsn);
            
            for gsn = 1:maxGsn
                GMmodels{gsn} = fitgmdist(projSpikes{el}(:,1:nDims),gsn,...
                    'Options',statset('MaxIter',500),...
                    'Start','plus','RegularizationValue',0.001);
                aic(gsn) = GMmodels{gsn}.AIC;
                bic(gsn) = GMmodels{gsn}.BIC;
            end
            
            thr = 0.1;
            gsnBest = find(bic < (1-thr)*bic(end)+thr*bic(1),1);
            
            clusterParams{el} = GMmodels{gsnBest};
            
            %% Assigning output
            neuronEls{el} = el*ones(gsnBest,1);
            neuronClusters{el} = (1:gsnBest)';
            
            spikeClust = clusterParams{el}.cluster(projSpikes{el}(:,1:nDims));
            
            spikeTimesNeuron{el} = cell(gsnBest,1);
            for gsn = 1:gsnBest
                spikeTimesNeuron{el}{gsn} = spikeTimesEl{el}(spikeClust == gsn);
            end
            
        catch error
            error
            disp(['Error at electrode ',num2str(el),', skipping.']);
        end
        
        if clustConfig.debug
            disp(sprintf(['Electrode ',num2str(el),':\n',...
            num2str(gsnBest),' neurons found.\n',...
            'Time for Gaussian Mixture Clustering ',num2str(toc(glmTimer)),' seconds\n',...
            '-----------------------\n']));
            
            if false % Debug plots
                %%
                figure(1)
                plotyy(1:maxGsn,bic,1:maxGsn,aic);
                legend('BIC','AIC');
                
                for gsn = 1:maxGsn
                    
                    GMmodel = GMmodels{gsn};
%                   idx = GMmodel.cluster(projSpikes{el}(:,1:nDims));
                    idx = (GMmodel.posterior(projSpikes{el}(:,1:nDims)) > 0.6)*(1:gsn)';
                    
                    
                    figure(2)
                    scatter3(projSpikes{el}(:,1),projSpikes{el}(:,2),projSpikes{el}(:,3),9,idx);
                    colormap jet
                    colorbar
                    title('Gaussian mixture')
                    
                    for g = 1:gsn
                        covMat = GMmodel.Sigma(1:3,1:3,g);
                        center = GMmodel.mu(g,1:3);
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
                    end
                    gsnBest
                    gsn
                end
            end % debug plots
        end % if debug
        
    end % el
    
    % Concatenating all neurons found
    neuronEls = vertcat(neuronEls{:});
    neuronClusters = vertcat(neuronClusters{:});
    spikeTimesNeuron = vertcat(spikeTimesNeuron{:});
    
end % function
