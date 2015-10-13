function [adjacent,maxAdjacent] = catchAdjWJava( dataSource, electrodeUsage )
    %CATCHADJWJAVA Takes a dataFileupsampler or an eiFile as input - outputs the adjacency list of electrodes
    %   Wraps the use of ElectrodeMapFactory class from other Matlab functions
    %   Output is numbered in Matlab Format: [1 - nElectrodes]
    % 
    % eiFile as input is a convenient patch for duplicate removal based on neighboring eis.
    % on a system without the .bin files to create a dataFileUpsampler.
    %
    % Author -- Vincent Deo -- Stanford University -- Oct 1st, 2015
    
    validateattributes(dataSource,{'DataFileUpsampler','edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile'},{},'','dataSource');
    
    import edu.ucsc.neurobiology.vision.electrodemap.*
    import edu.ucsc.neurobiology.vision.io.*
    
    % Get electrode map
    if isa(dataSource,'DataFileUpsampler');
        nElectrodes = dataSource.nElectrodes;
        header = dataSource.rawDataFile.getHeader();
        packedArrayID = int32(header.getArrayID());
    else
        nElectrodes = dataSource.nElectrodes;
        packedArrayID = int32(dataSource.getArrayID());
    end
    
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

