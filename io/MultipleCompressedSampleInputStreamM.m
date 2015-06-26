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
        % Two constructors implemented - Note there are 3 java constructors with added complexity
        function obj = MultipleCompressedSampleInputStreamM(varargin)
            % arguments (rawDataSource)
            % arguments (rawDataSource, bufferSizeInBytes, nBuffers)
            % arguments (rawDataSource, bufferSizeInBytes, nBuffers, waitForData)
            
            narginchk(1,3); % And should throw error further on if nargin = 2.
            
            import edu.ucsc.neurobiology.vision.io.*
            % TODO implement path validity check
            % rawDataSource = varargin{1} validity check
            if false
                ME = MException('MultipleCompressedSampleInputStreamM:IllegalConstructorArgumentException',...
                    'Constructor argument rawDataSource is not a valid data source.');
                throw(ME);
            end
            
            obj.rawDataSource = varargin{1};
            
            if nargin == 1
                obj.multipleCompressedSampleInputStream = ...
                    MultipleCompressedSampleInputStream(varargin{1});
            else
                validateattributes(varargin{2},{'numeric'},{'scalar','integer'},'','BufferSizeInBytes',2);
                validateattributes(varargin{3},{'numeric'},{'scalar','integer'},'','nBuffers',3);
                obj.multipleCompressedSampleInputStream = ...
                    MultipleCompressedSampleInputStream(varargin{1},varargin{2},varargin{3});
            end
            
        end % Constructor
        
        function header = getJavaHeader(obj)
            header = obj.multipleCompressedSampleInputStream.getHeader();
        end % getJavaHeader
        
        function finishSampleProcessing(obj)
            obj.multipleCompressedSampleInputStream.finishSampleProcessing();
        end % finishSampleProcessing
        
        function n = getCumulativeNSamples(obj)
            n = obj.multipleCompressedSampleInputStream.getCumulativeNSamples();
        end % getCumulativeNSamples
        
        function out = isFinished(obj)
           out = obj.multipleCompressedSampleInputStream.isFinished(); 
        end % getSamplesRead
        
        function n = getBufferUsage(obj)
            n = obj.multipleCompressedSampleInputStream.getBufferUsage();
        end % getNelectrodes
        
        function initialize(obj, buffer)
            % Need to map buffer into a java array for reference passing... buuuuuuut no
        end % initialize
        
        function processSample(obj, sample)
            % Need to map buffer into a java array for reference passing... buuuuuuut no
        end % processSample
        
        function start(obj)
            obj.multipleCompressedSampleInputStream.start();
        end % start
        
        function setInitializationListener(obj)
            % TODO Need to set up SampleInitializationListener interface for given matlab objects
        end
        
        function addSampleListener(obj)
            % TODO Need to set up SampleListener interface for given matlab objects
        end
        
        function removeSampleListener(obj)
            % TODO Need to set up SampleListener interface for given matlab objects
        end
        
        function n = getSamplesRead(obj)
            n = obj.multipleCompressedSampleInputStream.getSamplesRead();
        end
        
        function el = getNElectrodes(obj)
            el = obj.multipleCompressedSampleInputStream.getNElectrodes();
        end
        
    end % Methods
    
end % Class

