classdef folderHandler < handle
    %FOLDERHANDLER Class to handles files for the cluster visualizer backend
    % All constructor file checking, error handling to be done here
    % Files and file existence are private.
    % Backend constructor should proceed as
    %       folderHandler.hasFile(fh.FILE_TAG)
    %       folderHandler.getPath(fh.FILE_TAG)
    
    enumeration
        NEURON_MAT_FILE
        PRJ_MAT_FILE
        MODEL_MAT_FILE
        SPIKES_MAT_FILE
        CLEAN_MAT_FILE
        EI_FILE
        STA_FILE
        GLOBALS_FILE
        CLASSIFICATION_FILE
    end
    
    methods
    end
    
    methods (Static)
        function ext = getExtension(fileTag)
            switch fileTag
                case folderHandler.NEURON_MAT_FILE
                    ext = '.neurons.mat';
                case folderHandler.PRJ_MAT_FILE
                    ext = '.prj.mat';
                case folderHandler.MODEL_MAT_FILE
                    ext = '.model.mat';
                case folderHandler.SPIKES_MAT_FILE
                    ext = '.spikes.mat';
                case folderHandler.CLEAN_MAT_FILE
                    ext = '.clean.mat';
                case folderHandler.EI_FILE
                    ext = '.ei';
                case folderHandler.STA_FILE
                    ext = '.sta';
                case folderHandler.GLOBALS_FILE
                    ext = '.globals';
                case folderHandler.CLASSIFICATION_FILE
                    ext = '.txt';
                otherwise
                    throw(MException('','folderHandler:getExtension - invalid tag.'));
            end
        end
    end
    
end

