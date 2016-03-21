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
    fastRows = zeros(max(allIDs),1);
    [~,~,positions] = intersect(1:max(allIDs),allIDs);
    fastRows(allIDs) = positions;
    
    for i = 1:size(manualActions,1)
        action = manualActions{i,1};
        params = manualActions{i,2};
        data = manualActions{i,3};
        switch action
            case EditAction.CONSOLIDATE
                % Skip silently
            case EditAction.ELEVATE
                % find the rows
                r = fastRows(params{1});
                elevatedStatus(r) = true;
            case EditAction.MERGE
                % Hard merge
            case EditAction.SHRINK
                % Shrink. Reaffect spike train from data field.
            case EditAction.RECLUSTER
                % Recluster. Reaffect from data field ?
        otherwise
                throw(MException('','applyManualCleaningPattern - Unhandled EditAction in switch statement.'));
        end
    end
end

