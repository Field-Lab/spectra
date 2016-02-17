function concatenatedDuplicateRemoval( datasets, saveRoot, saveFolderAndName, timeCommands )
    %CONCATENATEDDUPLICATEREMOVAL Applies duplicate removal for a concatenated analysis
    %
    % This includes the following steps:
    %   Removing low-count and contaminated neurons from the concatenated neurons.mat
    %   Splitting the concatenated neurons.mat into the subfolders.
    %   Computing EIs for each sub-dataset
    %   merging the EIs using Join
    %
    % This function deals with all file management, and no returns are made to
    % ConcatenatedAnalysis (intended caller so far).
    %
    % Vincent Deo - Stanford University - Nov 9th 2015
    
    %% Low count and contam in concatenated .neurons.mat
    % Load concat - var neuronSpikeTimes, neuronEls, neuronClusters
    load([saveRoot,filesep,'concat.neurons.mat']);
    nNeurons = numel(neuronEls);
    %%%
    fprintf('Initial neurons: %u\n',nNeurons);
    %%%
    
    % get data length
    load([saveRoot,filesep,'concat.spikes.mat'],'nSamples');
    
    % Load Configuration
    cfg = mVisionConfig();
    cleanConfig = cfg.getCleanConfig();
    
    % Remove low count and contaminated neurons - augment low count by number of datasets factor
    toRemove = cellfun(@(x) numel(x) < cleanConfig.minSpikes * numel(saveFolderAndName),...
        neuronSpikeTimes);
    for i = 1:nNeurons
        if (toRemove(i) || edu.ucsc.neurobiology.vision.anf.NeuronCleaning.getContam(neuronSpikeTimes{i}, nSamples) > cleanConfig.maxCont)
            toRemove(i) = true;
        end
    end
    % Save IDs to clear - superfluous as concat.neurons.mat is overwritten below, nevermind
    IDsRemovedAtContam = NeuronSaverM.getIDs(neuronEls(toRemove),neuronClusters(toRemove));
    save([saveRoot,filesep,'concat.clean.mat'],'IDsRemovedAtContam');
    
    % Clear bad neurons
    neuronEls = neuronEls(~toRemove);
    neuronClusters = neuronClusters(~toRemove);
    neuronSpikeTimes = neuronSpikeTimes(~toRemove);
    nNeurons = size(neuronEls,1);
    toRemove = false(nNeurons,1);
    
    % Update concat.neurons.mat
    save([saveRoot,filesep,'concat.neurons.mat'],'neuronEls','neuronClusters','neuronSpikeTimes');
    
    %%%
    fprintf('After low count and high contamination: %u\n',nNeurons);
    %%%
    
    %% Split neurons file
    % splitNeurons behavior is to discard neurons that end up empty during the splitting
    splitNeurons(saveRoot, saveFolderAndName, timeCommands);
    
    %% Compute EIs for all sub-datasets
    % using linux process spawning + tcsh trick to resync
    command = '';
    sepchar = ':';
    cmgr = 'edu.ucsc.neurobiology.vision.calculations.CalculationManager';
    vcfg = sprintf('.%svision%sconfig.xml',filesep,filesep);
    calc = '"Electrophysiological Imaging Fast"';
    s = '';
    for d = 1:numel(saveFolderAndName)
        [saveFolder,datasetName,~] = fileparts(saveFolderAndName{d});
        
        % Build a neurons file
        neuronSaver = NeuronSaverM(datasets{d}, saveFolder, datasetName,'',0);
        load([saveFolderAndName{d},'.neurons.mat']);
        neuronSaver.pushAllNeurons(neuronEls, neuronClusters, neuronSpikeTimes);
        neuronSaver.close();
        
        s = [s,...
            sprintf('java -Xmx4g -Xss100m -cp ".%1$svision%1$s%2$s.%1$svision%1$sVision.jar" %3$s -c %4$s %5$s %6$s %7$s %8$f %9$u %10$u %11$u %12$u',...
            filesep,sepchar,cmgr,vcfg,calc,saveFolder,datasets{d},...
            cleanConfig.EITC, cleanConfig.EILP, cleanConfig.EIRP,...
            cleanConfig.EISp, cleanConfig.EInThreads),' &;'];
    end
    % Effective EI computation
    system([s,'wait;']);
    % All EIs are computed at the end of the child wait command
    % this puts tcsh syntax to good use for OS level parallel processing
    
    %% Join EIs for concatenated analysis
    % This should cope with missing IDs in some datasets but not all, ie neurons that are not
    % present in all sub datasets.
    s = '';
    for d = 1:numel(saveFolderAndName)
        [saveFolder,~,~] = fileparts(saveFolderAndName{d});
        s = [s,saveFolder,';'];
    end
    s = ['"',s,'"'];
    % EIMerger is a custom vision class, which is now compiled inside the jar.
    system(sprintf('java -Xmx4g -Xss100m -cp ".%1$svision%1$s%2$s.%1$svision%1$sVision.jar" edu.ucsc.neurobiology.vision.io.EIMerger %3$s %4$s',...
        filesep,sepchar,s,saveRoot));
    
    %% EI based duplicate removal on concatenated neurons
    load([saveRoot,filesep,'concat.neurons.mat']);
    
    % EI access setup
    [~,rootName,~] = fileparts(saveRoot);
    system(sprintf('mv %1$s%2$s%3$s.ei %1$s%2$sconcat.ei',saveRoot,filesep,rootName));
    eiPath = [saveRoot, filesep,'concat.ei'];
    eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(eiPath);
    
    
    % Misc for dup removal algorithm
    eiSize = cleanConfig.EILP + cleanConfig.EIRP + 1;
    minWindowSize = cleanConfig.globMinWin(2) - cleanConfig.globMinWin(1);
    samplingVal = cleanConfig.globMinWin(1):0.01:cleanConfig.globMinWin(2);
    
    nElectrodes = eiFile.nElectrodes;
    [adjacent,~] = catchAdjWJava( eiFile, 2 );
    
    % Saving merge pattern - col 1 neuron kept - col 2 neuron merged and discarded
    IDsMerged = zeros(0,2);
    
    % Single Electrode removal and merges
    for el = 2:nElectrodes
        elNeurInd = find(neuronEls == el);
        if numel(elNeurInd) <= 1
            continue;
        end
        
        trace = zeros(numel(elNeurInd), (eiSize - minWindowSize) * numel(adjacent{el}));
        elNeuronIDs = arrayfun(@(x) neuronSaver.getNeuronID(neuronEls(x),neuronClusters(x)), elNeurInd);
        
        % EI load
        for n = 1:numel(elNeurInd)
            fullEI = eiFile.getImage(elNeuronIDs(n));
            ei = squeeze(fullEI(1,adjacent{el},:))';
            % Realignment of minima
            traceTemp = zeros(size(ei,1)-minWindowSize, size(ei,2));
            for neighborEl = 1:size(ei,2)
                interpolant = griddedInterpolant(1:size(ei,1),ei(:,neighborEl),'spline');
                [~, offset] = min(interpolant(samplingVal));
                traceTemp(:,neighborEl) = interpolant((1:size(traceTemp,1)) + (offset - 1) / 100);
            end
            trace(n,:) = traceTemp(:)';
        end % n
        
        % Distance computation
        parts = returnL2MergeClasses(trace, cleanConfig.eiThrW);
        
        for cc = 1:numel(parts)
            if numel(parts{cc}) > 1
                spikeCounts = cellfun(@(spikeTrain) numel(spikeTrain), neuronSpikeTimes(elNeurInd(parts{cc})));
                [~,b] = max(spikeCounts);
                bestNeuronIndex = elNeurInd(parts{cc}(b));
                
                toRemove(elNeurInd(parts{cc})) = true;
                toRemove(bestNeuronIndex) = false;
                
                % Merging Pattern
                mergePairs = NeuronSaverM.getIDs(repmat([el el],numel(parts{cc}),1),...
                    neuronClusters([repmat(bestNeuronIndex,numel(parts{cc}),1),elNeurInd(parts{cc})]));
                mergePairs(b,:) = [];
                IDsMerged = [IDsMerged;mergePairs];
                
                neuronSpikeTimes{bestNeuronIndex} = sort(horzcat(neuronSpikeTimes{elNeurInd(parts{cc})}));
            end
        end % cc
    end % el
    
    % Save IDs to clear
    save([saveRoot,filesep,'concat.clean.mat'],'IDsMerged','-append');
    
    % Updating contents
    neuronEls = neuronEls(~toRemove);
    neuronClusters = neuronClusters(~toRemove);
    neuronSpikeTimes = neuronSpikeTimes(~toRemove);
    nNeurons = size(neuronEls,1);
    toRemove = false(nNeurons,1);
    
    %%%
    fprintf('After single el merge pass: %u\n',nNeurons);
    %%%
    
    % Dual Neighbor electrodes duplicate removal and discard
    for el = 2:nElectrodes
        for el2Index = 2:numel(adjacent{el})
            el2 = adjacent{el}(el2Index);
            
            elNeurInd = find(or(neuronEls == el, neuronEls == el2));
            intersectEls = intersect(adjacent{el},adjacent{el2});
            
            if numel(elNeurInd) <= 1
                continue;
            end
            
            trace = zeros(numel(elNeurInd), (eiSize - minWindowSize) * numel(intersectEls));
            elNeuronIDs = arrayfun(@(x) neuronSaver.getNeuronID(neuronEls(x),neuronClusters(x)), elNeurInd);
            
            % EI load
            for n = 1:numel(elNeurInd)
                fullEI = eiFile.getImage(elNeuronIDs(n));
                ei = squeeze(fullEI(1,intersectEls,:))';
                % Realignment of minima
                traceTemp = zeros(size(ei,1)-minWindowSize, size(ei,2));
                for neighborEl = 1:size(ei,2)
                    interpolant = griddedInterpolant(1:size(ei,1),ei(:,neighborEl),'spline');
                    [~, offset] = min(interpolant(samplingVal));
                    traceTemp(:,neighborEl) = interpolant((1:size(traceTemp,1)) + (offset - 1) / 100);
                end
                trace(n,:) = traceTemp(:)';
            end % n
            
            % Distance computation
            parts = returnL2MergeClasses(trace, cleanConfig.eiThrG);
            
            for cc = 1:numel(parts)
                if numel(parts{cc}) > 1
                    spikeCounts = cellfun(@(spikeTrain) numel(spikeTrain), neuronSpikeTimes(elNeurInd(parts{cc})));
                    [~,b] = max(spikeCounts);
                    bestNeuronIndex = elNeurInd(parts{cc}(b));
                    
                    wasRemoved = toRemove(bestNeuronIndex);
                    toRemove(elNeurInd(parts{cc})) = true;
                    toRemove(bestNeuronIndex) = wasRemoved;
                end
            end % cc
        end % el2
    end % el
    
    
    % Save IDs to clear
    IDsDuplicatesRemoved = NeuronSaverM.getIDs(neuronEls(toRemove),neuronClusters(toRemove));
    if numel(IDsDuplicatesRemoved) == 0
        IDsDuplicatesRemoved = zeros(0,2);
    end
    save([saveRoot,filesep,'concat.clean.mat'],'IDsDuplicatesRemoved','-append');
    
    neuronEls = neuronEls(~toRemove);
    neuronClusters = neuronClusters(~toRemove);
    neuronSpikeTimes = neuronSpikeTimes(~toRemove);
    nNeurons = size(neuronEls,1);
    %%%
    fprintf('Final neurons: %u\n',nNeurons);
    %%%
    
    % Update concat.neurons.mat
    save([saveRoot,filesep,'concat.neurons.mat'],'neuronEls','neuronClusters','neuronSpikeTimes');
    
    %% Apply cleaning pattern to all sub datasets and convert to .neurons
    
    for d = 1:numel(datasets)
        load([saveFolderAndName{d},'.neurons.mat']);
        % Apply cleaning pattern
        [neuronEls, neuronClusters, neuronSpikeTimes] = ...
            applyCleaningPattern( IDsRemovedAtContam, IDsMerged, IDsDuplicatesRemoved,...
            neuronEls, neuronClusters, neuronSpikeTimes );
        save([saveFolderAndName{d},'.neurons.mat'],...
            'neuronEls','neuronClusters','neuronSpikeTimes');
        
        % Matlab to vision neuron file convert
        [saveFolder,datasetName,~] = fileparts(saveFolderAndName{d});
        saver = NeuronSaverM(datasets{d},saveFolder,datasetName,'',0);
        saver.pushAllNeurons(neuronEls, neuronClusters, neuronSpikeTimes);
        saver.close();
    end
end

