classdef MultipleCompressedSampleInputStreamM
    %MULTIPLECOMPRESSEDSAMPLEINPUTSTREAMM Wrapper of the java class MultipleCompressedInputStream
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable, GetAccess = protected)
        multipleCompressedSampleInputStream
    end
    properties (SetAccess = immutable, GetAccess = public)
        rawDataSource@char % Output folder in String format - provided to constructor and checked
    end
    
    methods
        % Basic constructor0 implemented - Note there are 3 java constructors with added complexity
        function obj = MultipleCompressedSampleInputStreamM(rawDataSource)
            import edu.ucsc.neurobiology.vision.io.*
            % TODO implement path validity check
            % rawDataSource validity check
            if false 
               ME = MException('MultipleCompressedSampleInputStreamM:IllegalConstructorArgumentException',...
                    'Constructor argument rawDataSource is not a valid data source.');
                throw(ME); 
            end
            obj.rawDataSource = rawDataSource;
            % Call java constructor and store
            obj.multipleCompressedSampleInputStream = MultipleCompressedSampleInputStream(rawDataSource);          
        end        
        
    end
    
end

