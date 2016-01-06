function CC = returnL2MergeClasses( arrays, threshold )
    %BUILPAIRWISEL2 Computes matrix of (almost) pairwise L2 distance between a set of
    % vectors - then processes the list of connected components for a given splitting threshold
    %
    % Metric on which thresholding is made: some sort of chi2
    %    || a - b ||^2
    % ------------------
    % || a || .* || b ||
    %
    % arrays - row oriented input
    % threshold - threshold to use
    %
    % CC - connected compents in a column cell array - each cell being a row of all indices to
    % merge
    
    % Distance matrix
    norms = sqrt(sum(arrays.^2,2));
    dists = sum(bsxfun(@minus,permute(arrays,[3,1,2]),permute(arrays,[1,3,2])).^2,3);
    dists = bsxfun(@rdivide,bsxfun(@rdivide,dists,norms),norms');
    
    % Process tree
    dists(isnan(dists)) = 0;
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
        g = Graph(max(max(pairList(:,1:2))));
    
    % Process all edges
    g.addAllEdges(pairList(:,1)-1,pairList(:,2)-1,pairList(:,3));
    
    % Extract final connected component structure and merge edge values list
    % g.removeSingletons();
    % t = g.getFinalTree();
    % CC = t.partitioning(threshold);
    CC = g.getFinalForest();
    
    if ~isa(CC,'cell') % java return may be array type if it has appropriate dimensions
        CC = mat2cell(CC,ones(1,size(CC,1)));
    end
    
    CC = cellfun(@(x) x+1,CC,'uni',false);
end

