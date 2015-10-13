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
    % Author -- Vincent Deo -- Stanford University -- October 12, 2015
    
    
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
    
    % Remove low count and contaminated neurons
    toRemove = cellfun(@(x) numel(x) < cleanConfig.minSpikes, neuronSpikeTimes);
    for i = 1:nNeurons
        if (toRemove(i) || edu.ucsc.neurobiology.vision.anf.NeuronCleaning.getContam(neuronSpikeTimes{i}, nSamples) > cleanConfig.maxCont)
            toRemove(i) = true;
        end
    end
    
    % Clear bad neurons
    neuronEls = neuronEls(~toRemove);
    neuronClusters = neuronClusters(~toRemove);
    neuronSpikeTimes = neuronSpikeTimes(~toRemove);
    nNeurons = size(neuronEls,1);
    toRemove = false(nNeurons,1);
    
    
    %%%
    fprintf('After low count and high contamination: %u\n',nNeurons);
    %%%
    
    
    % Build a neurons file and compute EIs
    neuronSaver = NeuronSaverM(dataPath, saveFolder, datasetName,'');
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
    
    system(sprintf('java -cp ".%1$svision%1$s%2$s.%1$svision%1$sVision.jar" %3$s -c %4$s %5$s %6$s %7$s %8$f %9$u %10$u %11$u %12$u',...
        filesep,sepchar,cmgr,vcfg,calc,saveFolder,dataPath,...
        cleanConfig.EITC, cleanConfig.EILP, cleanConfig.EIRP,...
        cleanConfig.EISp, cleanConfig.EInThreads));
    delete([saveFolder, filesep, datasetName,'.neurons']);

    % EI access setup
    eiPath = [saveFolder, filesep, datasetName,'.ei'];
    eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(eiPath);
    
    
    % Misc for dup removal algorithm
    eiSize = cleanConfig.EILP + cleanConfig.EIRP + 1;
    minWindowSize = cleanConfig.globMinWin(2) - cleanConfig.globMinWin(1);
    samplingVal = cleanConfig.globMinWin(1):0.01:cleanConfig.globMinWin(2);
    
    nElectrodes = eiFile.nElectrodes;
    [adjacent,~] = catchAdjWJava( eiFile, 2 );
    
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
                
                neuronSpikeTimes{bestNeuronIndex} = sort(horzcat(neuronSpikeTimes{elNeurInd(parts{cc})}));
            end
        end % cc
    end % el
    
    
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
            parts = returnL2MergeClasses(trace, cleanConfig.eiThrA);
            
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
    
    neuronEls = neuronEls(~toRemove);
    neuronClusters = neuronClusters(~toRemove);
    neuronSpikeTimes = neuronSpikeTimes(~toRemove);
    nNeurons = size(neuronEls,1);
    %%%
    fprintf('Final neurons: %u\n',nNeurons);
    %%%
end

