function trackPerm = optimalBlockDiagPerm( matrix )
    %OPTIMALBLOCKDIAGPERM Returns optimal permutation to give a block diagonal
    %aspect to a symmetric matrix.
    %
    % Proceeds by greedily finding transpositions reducing a weighted cost
    % strongly penalizing values far from the diagonal
    
    n = size(matrix,1);
    m = matrix ./ max(abs(max(max(matrix))),abs(min(min(matrix))));
    
    weight = toeplitz((0:(n-1)).^1);
    currValue = sum(sum(weight.*abs(m)));
    trackPerm = 1:n;
    for k = 1:(2*n)
        transpoValue = Inf(n,n);
        for i = 1:(n-1)
            for j = i:n
                perm = 1:n;
                perm(i) = j;
                perm(j) = i;
                transpoValue(i,j) = sum(sum(weight.*m(perm,perm).^2));
            end
        end
        [i,j] = find(transpoValue == min(min(transpoValue)));
        perm = 1:n;
        perm(i) = j;
        perm(j) = i;
        m = m(perm,perm);
        trackPerm = trackPerm(perm);
        currValue = transpoValue(i,j);
    end
end

