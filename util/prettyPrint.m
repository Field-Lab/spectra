function s = prettyPrint(vect)
    % function prettyPrint
    %   Generates a space separated, bracket delimited string for an integer array.

    if size(vect,2) == 1
        vect =vect';
    end
    s = ['[ ',num2str(vect,'%u '),']'];
end

