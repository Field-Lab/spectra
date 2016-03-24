function s = prettyPrint(vect)
    % function prettyPrint
    %   Generates a space separated, bracket delimited string for an integer array.

    if size(vect,2) == 1
        vect =vect';
    end
    if numel(vect) == 0
        s = '[]';
        return;
    end
    
    if nnz(vect - round(vect)) == 0 % integer
        s = ['[ ',num2str(vect,'%u '),' ]'];
    else
        s = ['[ ',num2str(vect,'%.2f '),' ]'];
    end
end

