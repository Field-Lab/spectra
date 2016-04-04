function [neuronEls, neuronClusters, neuronSpikeTimes, elevatedStatus] = ...
        applyManualCleaningPattern(neuronEls, neuronClusters,...
        neuronSpikeTimes, manualActions, elevatedStatus)
    
    nNeurons = size(neuronEls,1);
    
    % Argument check
    validateattributes(neuronEls,{'numeric'},{'size',[nNeurons, 1]},'','neuronEls');
    validateattributes(neuronClusters,{'numeric'},{'size',[nNeurons, 1]},'','neuronClusters');
    validateattributes(neuronSpikeTimes,{'cell'},{'size',[nNeurons, 1]},'','neuronSpikeTimes');
    validateattributes(manualActions,{'cell'},{'size',[nan, 3]},'','manualActions');
    validateattributes(elevatedStatus,{'logical'},{'size',[nNeurons, 1]},'','elevatedStatus');
    
    allIDs = NeuronSaverM.getIDs(neuronEls, neuronClusters);
    fastRows = zeros(max(allIDs) + 15,1);
    % " + 15 ": new IDs may be added on the last electrode. Conservative.
    [~,~,positions] = intersect(1:max(allIDs),allIDs);
    fastRows(allIDs) = positions;
    
    % removal flags
    toRemove = false(nNeurons,1);
    
    for i = 1:size(manualActions,1)
        action = manualActions{i,1};
        params = manualActions{i,2};
        data = manualActions{i,3};
        
        % Manual actions should come in order
        % Throw errors rather than recover in case of missing IDs
        % User decides: we do not care about elevated neurons
        switch action
            case EditAction.CONSOLIDATE
                % Skip silently
            case EditAction.KEEP
                % find the rows
                r = fastRows(params{1});
                elevatedStatus(r) = true;
            case EditAction.DELETE
                r = fastRows(params{1});
                elevatedStatus(r) = true;
                toRemove(r) = true;
            case EditAction.MERGE
                % Hard merge
                % find the rows
                r = fastRows(params{1});
                if any(r == 0)
                    throw(MException('','applyManualCleaningPattern - IDs requested for merge are missing'));
                end
                % Find master ID
                nSpikes = cellfun(@(x) numel(x),neuronSpikeTimes(r),'uni',true);
                [~,b] = max(nSpikes);
                rMaster = r(b);
                r(b) = [];
                toRemove(r) = true;
                neuronSpikeTimes{rMaster} = sort(horzcat(neuronSpikeTimes{[rMaster;r]}),'ascend');
            case EditAction.SHRINK
                r = fastRows(params{1});
                neuronSpikeTimes(r) = data{3}; % Forwarded spike trains
                elevatedStatus(r) = true;
            case EditAction.RECLUSTER
                r = fastRows(params{1});
                toRemove(r) = true;
                newN = numel(data{1});
                [el, clust] = NeuronSaverM.getElClust(data{1});
                % Append, and update fastrows
                fastRows(data{1}(:)) = numel(allIDs) + (1:newN);
                allIDs = [allIDs ; data{1}(:)];
                neuronEls = [neuronEls; el(:)];
                neuronClusters = [neuronClusters; clust(:)];
                neuronSpikeTimes = [neuronSpikeTimes ; data{2}(:)];
                toRemove = [toRemove ; false(newN, 1)];
                elevatedStatus = [elevatedStatus ; true(newN, 1)];
        otherwise
                throw(MException('','applyManualCleaningPattern - Unhandled EditAction in switch statement.'));
        end
    end
    
    % Remove discarded rows
    allIDs(toRemove) = [];
    neuronEls(toRemove) = [];
    neuronClusters(toRemove) = [];
    neuronSpikeTimes(toRemove) = [];
    elevatedStatus(toRemove) = [];
    
    % Re-sort
    [~, idx] = sort(allIDs);
    neuronEls = neuronEls(idx);
    neuronClusters = neuronClusters(idx);
    neuronSpikeTimes = neuronSpikeTimes(idx);
    elevatedStatus = elevatedStatus(idx);

end

