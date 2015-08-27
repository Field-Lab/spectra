function [adjacent,maxAdjacent] = catchAdjWJava( dataSource, electrodeUsage )
    %CATCHADJWJAVA Takes a dataFileupsampler as input - outputs the adjacency list of electrodes
    %   Wraps the use of ElectrodeMapFactory class from other Matlab functions
    %   Output is numbered in Matlab Format: 1 - nElectrodes
    % 
    % Author -- Vincent Deo -- Stanford University -- August 27, 2015
    
    validateattributes(dataSource,{'DataFileUpsampler'},{},'','dataSource');
    
    import edu.ucsc.neurobiology.vision.electrodemap.*
    import edu.ucsc.neurobiology.vision.io.*
    
    % Get electrode map
    nElectrodes = dataSource.nElectrodes;
    header = dataSource.rawDataFile.getHeader();
    packedArrayID = int32(header.getArrayID());
    electrodeMap = ElectrodeMapFactory.getElectrodeMap(packedArrayID);
    
    adjacent = cell(nElectrodes,1);
    maxAdjacent = 0;
    
    % Get adjacents
    for el = 0:(nElectrodes-1)
        adjacent{el+1} = electrodeMap.getAdjacentsTo(el, electrodeUsage) + 1 ;
        if numel(adjacent{el+1}) > maxAdjacent
            maxAdjacent = numel(adjacent{el+1});
        end
    end % el
    
end

