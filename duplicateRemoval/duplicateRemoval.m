function [neuronEls, neuronClusters, neuronSpikeTimes] = ...
        duplicateRemoval( dataPath, saveFolder, datasetName, timeTag, neuronEls, neuronClusters, neuronSpikeTimes)
    % DUPLICATE REMOVAL - Remove duplicates neurons based on EI comparisons
    %
    % Removal of low count and above threshold contamination neurons
    % Removal within electrode of neurons with similar EIs
    % Removal across neighboring electrodes of neurons with similar EIs
    % EI comparison is based on Chi2-like distance metric on realigned timecourses
    % over relevant electrodes
    %
    % Inputs:
    %   dataPath: path to the data folder
    %   saveFolder: path to the analysis folder
    %   timeTag: additional time command "(***-***)" on which this analysis is being done
    %   neuronEls: nNeurons x 1 array containing the electrode of each neurons
    %   neuronClusters: nNeurons x 1 array containing the cluster number of each neuron
    %   neuronSpikeTimes: nNeurons x 1 cell array, with each cell a vector containing
    %       all the spike times for this neuron
    %
    % Return:
    %   neuronEls: same spec as input
    %   neuronClusters: same spec as input
    %   neuronSpikeTimes: same spec as input
    %
    % Saves:
    %   IDs of removed neurons and merged neurons, in file cleanPattern.mat
    %   For purpose of broadcasting the cleaning pattern in concatenated analyses.
    %
    % Author -- Vincent Deo -- Stanford University -- October 12, 2015
    %% Setup
    
    nNeurons = size(neuronEls,1);
    validateattributes(neuronEls,{'numeric'},{'size',[nNeurons, 1]},'','neuronEls');
    validateattributes(neuronClusters,{'numeric'},{'size',[nNeurons, 1]},'','neuronClusters');
    validateattributes(neuronSpikeTimes,{'cell'},{'size',[nNeurons, 1]},'','neuronSpikeTimes');
    
    %%%
    fprintf('Initial neurons: %u\n',nNeurons);
    %%%
    
    % get data length
    datasource = DataFileUpsampler([dataPath,timeTag]);
    nSamples = datasource.stopSample - datasource.startSample; % horribly ugly - won't go through concatenation patch, etc...
    datasource = [];
    
    % Load Configuration
    cfg = mVisionConfig();
    cleanConfig = cfg.getCleanConfig();
    
    %% Remove low count and contaminated neurons
    toRemove = cellfun(@(x) numel(x) < cleanConfig.minSpikes, neuronSpikeTimes);
    for i = 1:nNeurons
        if (toRemove(i) ||...
                edu.ucsc.neurobiology.vision.anf.NeuronCleaning.getContam(neuronSpikeTimes{i}, nSamples) >...
                cleanConfig.maxCont)
            toRemove(i) = true;
        end
    end
    
    % Save IDs to clear
    IDsRemovedAtContam = NeuronSaverM.getIDs(neuronEls(toRemove),neuronClusters(toRemove));
    save([saveFolder,filesep,'cleanPattern.mat'],'IDsRemovedAtContam');
    
    % Clear bad neurons
    neuronEls = neuronEls(~toRemove);
    neuronClusters = neuronClusters(~toRemove);
    neuronSpikeTimes = neuronSpikeTimes(~toRemove);
    nNeurons = size(neuronEls,1);
    
    % Reinitialize
    toRemove = false(nNeurons,1);
    
    %%%
    fprintf('After low count and high contamination: %u\n',nNeurons);
    %%%
    
    
    %% Build a neurons file and compute EIs
    neuronSaver = NeuronSaverM(dataPath, saveFolder, datasetName,'',0);
    neuronSaver.pushAllNeurons(neuronEls, neuronClusters, neuronSpikeTimes);
    neuronSaver.close();
    
    % EI computation
    if isunix
        sepchar = ':';
    else
        sepchar = ';';
    end
    cmgr = 'edu.ucsc.neurobiology.vision.calculations.CalculationManager';
    vcfg = sprintf('.%svision%sconfig.xml',filesep,filesep);
    calc = '"Electrophysiological Imaging Fast"';
    % Grind-style call for EIs computation at system level
    system(sprintf('java -Xmx4g -Xss100m -cp ".%1$svision%1$s%2$s.%1$svision%1$sVision.jar" %3$s -c %4$s %5$s %6$s %7$s %8$f %9$u %10$u %11$u %12$u',...
        filesep,sepchar,cmgr,vcfg,calc,saveFolder,dataPath,...
        cleanConfig.EITC, cleanConfig.EILP, cleanConfig.EIRP,...
        cleanConfig.EISp, cleanConfig.EInThreads));
    delete([saveFolder, filesep, datasetName,'.neurons']);
    
    %% EI access setup
    eiPath = [saveFolder, filesep, datasetName,'.ei'];
    eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(eiPath);
    
    % Misc for dup removal algorithm
    eiSize = cleanConfig.EILP + cleanConfig.EIRP + 1;
    minWindowSize = cleanConfig.globMinWin(2) - cleanConfig.globMinWin(1);
    samplingVal = cleanConfig.globMinWin(1):cleanConfig.resamplePitch:cleanConfig.globMinWin(2);
    basicValues = 1:(eiSize - minWindowSize);
    
    nElectrodes = eiFile.nElectrodes;
    [adjacent1,~] = catchAdjWJava( eiFile, 1 );
    [adjacent2,~] = catchAdjWJava( eiFile, 2 );
    
    % EI loader
    eiStorage = cell(nNeurons,1); % full remaining neurons EIs
    isRealigned = false(nNeurons,nElectrodes); % tag to upsample any ID/electrode only once
    
    % Preload all EIs
    allIDs = NeuronSaverM.getIDs(neuronEls,neuronClusters);
    for n = 1:nNeurons
        eiStorage{n} = eiFile.getImage(allIDs(n));
        eiStorage{n} = squeeze(eiStorage{n}(1,:,:))';
    end
    fprintf('EIs loaded.\n');
    
    %% Single Electrode removal and merges
    % Saving merge pattern - col 1 neuron kept - col 2 neuron merged and discarded
    IDsMerged = [];
    
    for el = 2:nElectrodes
        elNeurInd = find(neuronEls == el);
        if numel(elNeurInd) <= 1
            continue;
        end
        
        trace = zeros(numel(elNeurInd), (eiSize - minWindowSize) * numel(adjacent2{el}));
        
        % EI upsampling handling
        traceTemp = zeros(eiSize - minWindowSize, numel(adjacent2{el}));
        for n = 1:numel(elNeurInd)
            % Realignment procedure
            for neighborElInd = 1:numel(adjacent2{el});
                neighborEl = adjacent2{el}(neighborElInd);
                if ~isRealigned(elNeurInd(n),neighborEl)
                    interpolant = griddedInterpolant(1:eiSize,eiStorage{elNeurInd(n)}(:,neighborEl),'spline');
                    [~, offset] = min(interpolant(samplingVal));
                    eiStorage{elNeurInd(n)}(basicValues,neighborEl) =...
                        interpolant(basicValues + (offset - 1) * cleanConfig.resamplePitch);
                    isRealigned(elNeurInd(n),neighborEl) = true;
                end
                traceTemp(:,neighborElInd) = eiStorage{elNeurInd(n)}(basicValues,neighborEl);
            end
            trace(n,:) = traceTemp(:)';
        end % n
        % L2 normalization
        trace = bsxfun(@rdivide,trace,sqrt(sum(trace.^2,2)));
        
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
    save([saveFolder,filesep,'cleanPattern.mat'],'IDsMerged','-append');
    
    % Updating contents
    neuronEls = neuronEls(~toRemove);
    neuronClusters = neuronClusters(~toRemove);
    neuronSpikeTimes = neuronSpikeTimes(~toRemove);
    nNeurons = size(neuronEls,1);
    
    allIDs = allIDs(~toRemove);
    eiStorage = eiStorage(~toRemove);
    isRealigned = isRealigned(~toRemove,:);
    
    % Reinitialize
    toRemove = false(nNeurons,1);
    
    %%%
    fprintf('After single electrode pair merges: %u\n',nNeurons);
    %%%
    
    %% Global duplicate removal and discard (inc. axonal detections of strong cells)
    % As we may need ALL realigned EI pieces, EI access strategy is changed
    for el = 2:nElectrodes
        for el2 = (el+1):nElectrodes
            elNeurInd1 = find(neuronEls == el);
            elNeurInd2 = find(neuronEls == el2);
            intersectEls = union(adjacent1{el},adjacent1{el2});
            
            if numel(elNeurInd1) == 0 || numel(elNeurInd2) == 0
                continue;
            end
            elNeurInd = [elNeurInd1;elNeurInd2];
            
            trace = zeros(numel(elNeurInd), (eiSize - minWindowSize) * numel(intersectEls));
            
            % EI upsampling handling
            traceTemp = zeros(eiSize - minWindowSize, numel(intersectEls));
            for n = 1:numel(elNeurInd)
                % Realignment procedure
                for neighborElInd = 1:numel(intersectEls);
                    neighborEl = intersectEls(neighborElInd);
                    if ~isRealigned(elNeurInd(n),neighborEl)
                        interpolant = griddedInterpolant(1:eiSize,eiStorage{elNeurInd(n)}(:,neighborEl),'spline');
                        [~, offset] = min(interpolant(samplingVal));
                        eiStorage{elNeurInd(n)}(basicValues,neighborEl) =...
                            interpolant(basicValues + (offset - 1) * cleanConfig.resamplePitch);
                        isRealigned(elNeurInd(n),neighborEl) = true;
                    end
                    traceTemp(:,neighborElInd) = eiStorage{elNeurInd(n)}(basicValues,neighborEl);
                end
                trace(n,:) = traceTemp(:)';
            end % n
            % L2 normalization
            trace = bsxfun(@rdivide,trace,sqrt(sum(trace.^2,2)));
            
            % Distance computation
            [parts,~] = returnL2MergeClasses(trace, cleanConfig.eiThrG);
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
    
    
    %% Save IDs to clear
    IDsDuplicatesRemoved = allIDs(toRemove);
    save([saveFolder,filesep,'cleanPattern.mat'],'IDsDuplicatesRemoved','-append');
    
    neuronEls = neuronEls(~toRemove);
    neuronClusters = neuronClusters(~toRemove);
    neuronSpikeTimes = neuronSpikeTimes(~toRemove);
    nNeurons = size(neuronEls,1);
    
    %%%
    fprintf('After neighbor electrode pair discard: %u\n',nNeurons);
    %%%
end

