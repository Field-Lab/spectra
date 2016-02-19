function parforGlobalizer( local_copy, name )
    %PARFORGLOBALIZER Summary of this function goes here
    %   Globalizes a matlab workspace global object
    %   To a parfor worker workspace
    %   Using a local copy to enter the parfor loop
    %   Then globalizing in this file, under name provided as second argument
    %   A deep copy is made when entering the parfor loop
    %
    % USE:
    %
    % global name
    % %% name is available
    % local_copy = name
    % parfor
    % %% name is unavailable
    % %% local_copy is still available
    % parforGlobalizer(local_copy,'name');
    % %% name is still unavailable
    % %% but a subfunction may now call 'global name'
    % %% and access the worker's workspace global copy.
    %
    % Author - Vincent Deo - Stanford University - Feb 18, 2016
    
    eval(sprintf('global %s;',name));
    eval(sprintf('%s = local_copy;',name));
end

