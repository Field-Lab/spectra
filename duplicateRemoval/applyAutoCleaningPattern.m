function [neuronEls, neuronClusters, neuronSpikeTimes] = ...
        applyAutoCleaningPattern(neuronEls, neuronClusters,...
        neuronSpikeTimes, autoActions, elevatedStatus)
    
    nNeurons = size(neuronEls,1);
    
    % Argument check
    validateattributes(neuronEls,{'numeric'},{'size',[nNeurons, 1]},'','neuronEls');
    validateattributes(neuronClusters,{'numeric'},{'size',[nNeurons, 1]},'','neuronClusters');
    validateattributes(neuronSpikeTimes,{'cell'},{'size',[nNeurons, 1]},'','neuronSpikeTimes');
    validateattributes(autoActions,{'cell'},{'size',[nan, 3]},'','autoActions');
    validateattributes(elevatedStatus,{'logical'},{'size',[nNeurons, 1]},'','elevatedStatus');
    listElevatedRows = find(elevatedStatus);
    
    allIDs = NeuronSaverM.getIDs(neuronEls, neuronClusters);
    fastRows = zeros(max(allIDs) + 15,1);
    [~,~,positions] = intersect(1:max(allIDs),allIDs);
    fastRows(allIDs) = positions;
    
    % removal flags
    toRemove = false(nNeurons,1);
    
    % Apply all automatic actions
    for i = 1:size(autoActions,1)
        action = autoActions{i,1};
        params = autoActions{i,2};
        
        switch action
            case EditAction.AUTO_RM
                % find the rows
                r = fastRows(params{1}); r(r == 0) = [];
                % remove the elevated rows
                r = setdiff(r,listElevatedRows);
                % discard
                toRemove(r) = true;
            case EditAction.AUTO_MERGE
                % find the rows
                rMaster = fastRows(params{1});
                r = fastRows(params{2}); r(r == 0) = [];
                % remove elevated rows
                r = setdiff(r,listElevatedRows);
                % depends on whether master ID is elevated
                if ~(rMaster == 0) && ~elevatedStatus(rMaster)
                    toRemove(r) = true;
                    neuronSpikeTimes{rMaster} = sort(horzcat(neuronSpikeTimes{[rMaster;r]}),'ascend');
                else % master is elevated, need to find another master neuron
                    if numel(r) > 0
                        nSpikes = cellfun(@(x) numel(x),neuronSpikeTimes(r),'uni',true);
                        [~,b] = max(nSpikes);
                        rMaster = r(b);
                        r(b) = [];
                        toRemove(r) = true;
                        neuronSpikeTimes{rMaster} = sort(horzcat(neuronSpikeTimes{[rMaster;r]}),'ascend');
                    end
                end
            case EditAction.AUTO_RM_DUP
                % find the rows
                r = fastRows(params{2});  r(r == 0) = [];
                % remove elevated rows
                r = setdiff(r,listElevatedRows);
                % depends on whether master ID is elevated
                toRemove(r) = true;
            otherwise
                throw(MException('','applyAutoCleaningPattern - Unhandled EditAction in switch statement.'));
        end
    end
    
    % Remove discarded rows
    neuronEls(toRemove) = [];
    neuronClusters(toRemove) = [];
    neuronSpikeTimes(toRemove) = [];
end

