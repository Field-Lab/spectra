function [CC,varargout] = returnL2MergeClasses( arrays, threshold )
    %BUILPAIRWISEL2 Computes matrix of (almost) pairwise L2 distance between a set of
    % vectors - then processes the list of connected components for a given splitting threshold
    %
    % Metric on which thresholding is made: L_2 norm
    %    || a - b ||^2
    % Please pre normalize input if using normalized L_2, using eg:
    % arrays = bsxfun(@rdivide,arrays,sqrt(sum(arrays.^2,2)));
    %
    % arrays - row oriented input
    % threshold - threshold to use
    %
    % CC - connected compents in a column cell array - each cell being a row of all indices to
    % merge
    
    % Distance matrix
    dists = sum(bsxfun(@minus,permute(arrays,[3,1,2]),permute(arrays,[1,3,2])).^2,3);
    
    % Process tree
    dists(isnan(dists)) = 0;
    varargout{1} = dists;
    [row,col] = find( and(dists > 0,dists < threshold) ); % remove singularities and values above threshold
    % remove diagonal terms if any
    diagRem = row == col; row(diagRem) = []; col(diagRem) = [];
    % remove duplicate indices for computing each edge only once
    duplicatePairs = unique(sort([row,col],2),'rows');
    
    % List of distance-valued connected neuron pairs
    % [row, col, edge_value]
    pairList = [double(duplicatePairs),dists(sub2ind(size(dists),duplicatePairs(:,1),duplicatePairs(:,2)))];
    % Sort by increasing edge value - required for graph/conn component analysis
    pairList = sortrows(pairList,3);

    % Threshold-based connected component analysis
    nodes = size(arrays,1);
    g = Graph(nodes);
    
    % Process all edges
    if numel(pairList) > 0
        g.addAllEdges(pairList(:,1)-1,pairList(:,2)-1,pairList(:,3));
    end
    % Extract final connected component structure and merge edge values list
    CC = g.getFinalForestEmptyRow();
    
    % Checking if C is cell type and using mat2cell if not is very slow.
    % getFinalForestEmptyRow() returns the final forest with added null row
    % to force matlab to take as cell array type (row dimension mismatch).
    % We discard this extra row:
    CC = CC(1:(end-1));
    
    CC = cellfun(@(x) x+1,CC,'uni',false);
    
end

