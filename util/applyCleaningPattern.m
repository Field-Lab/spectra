function [neuronEls, neuronClusters, neuronSpikeTimes] = ...
        applyCleaningPattern( IDsRemovedAtContam, IDsMerged, IDsDuplicatesRemoved,...
        neuronEls, neuronClusters, neuronSpikeTimes )
    %APPLYCLEANINGPATTERN Applies the neuron deduplication pattern made on 1 file
    % to another neuron file.
    %
    % Allows to broadcast a deduplication pattern that was found for the EIs of
    % a single file to the other files of a concatenated run.
    %
    % Inputs:
    %       IDsRemovedAtContam: neuron IDs of all the neurons discarded by the contamination
    %           threshold in the reference run.
    %       IDsMerged: IDs of neurons that were merged at the single electrode
    %           EI pass on the reference run.
    %       IDsDuplicatesRemoved: neuron IDs that were removed during the 19
    %           electrode EI pass on the reference run
    %   Note: the three arrays abov were processed in this order in the reference run
    %       but are supposed to be exclusive, and can be applied in any order.
    %       neuronEls: electrodes of neurons to clean
    %       neuronClusters: cluster numbers of neurons to clean
    %       neuronSpikeTrains: spike trains of neurons to clean
    %
    % Returns:
    %       neuronEls: cleaned neuron electrodes
    %       neuronClusters: cleaned neuron clusterNumbers
    %       neuronSpikeTimes: cleaned neurons spikeTrains
    %
    % Author -- Vincent Deo -- Stanford University -- November 2, 2015
    
    nNeurons = size(neuronEls,1);
    
    % Argument check
    validateattributes(neuronEls,{'numeric'},{'size',[nNeurons, 1]},'','neuronEls');
    validateattributes(neuronClusters,{'numeric'},{'size',[nNeurons, 1]},'','neuronClusters');
    validateattributes(neuronSpikeTimes,{'cell'},{'size',[nNeurons, 1]},'','neuronSpikeTimes');   
    
    allIDs = NeuronSaverM.getIDs(neuronEls, neuronClusters);
    
    % Find merge indexes - error will happen if neurons IDs we're seeking are missing
    [~,~,mergeKeep] = intersect(IDsMerged(:,1),allIDs);
    [~,~,mergeDiscard] = intersect(IDsMerged(:,2),allIDs);
    
    % Apply merges
    for i = 1:size(mergeKeep,1)
        neuronSpikeTimes{mergeKeep(i)} = union(neuronSpikeTimes{mergeKeep(i)},...
            neuronSpikeTimes{mergeDiscard(i)},'sorted');
    end
    
    % Remove non-merge neurons
    if size(IDsDuplicatesRemoved,2) == 2 % Stored removal tracking
        [~,~,ind] = intersect(union(union(IDsRemovedAtContam,IDsDuplicatesRemoved(:,2)),IDsMerged(:,2)), allIDs);    
    else % Did not store removal tracking
        [~,~,ind] = intersect(union(union(IDsRemovedAtContam,IDsDuplicatesRemoved(:,1)),IDsMerged(:,2)), allIDs);
    end
    neuronEls(ind) = [];
    neuronClusters(ind) = [];
    neuronSpikeTimes(ind) = [];
    
end

