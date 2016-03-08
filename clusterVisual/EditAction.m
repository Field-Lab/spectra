classdef EditAction
    %EDITACTION enumeration class describing possible edition actions
    %   This class specifies the possible edition action and
    %   Their required parameters are defined by the method checkParameterStructure
    %
    %   See ClusterEditBackend for how actions are applied
    %   and ClusterEditGUI for how they are called
    
    enumeration
        DEBUG
        DEBUG_2
    end
    
    methods
        
        function [valid,msg] = checkParameterStructure(obj,params)
            valid = false;
            msg = '';
            try
                switch obj
                    case EditAction.DEBUG
                        validateattributes(params,{'cell'},{'size',[1 1]},'','parameter cell array');
                        validateattributes(params{1},{'numeric'},{'row','nonempty'},'','DEBUG random nonempty list',1);
                    case EditAction.DEBUG_2
                        validateattributes(params,{'cell'},{'size',[2 1]},'','parameter cell array');
                        validateattributes(params{1},{'numeric'},{'row','nonempty'},'','DEBUG_2 random nonempty list',1);
                        validateattributes(params{2},{'char'},{'row','nonempty'},'','DEBUG_2 nonempty string',2);
                end
                valid = true;
            catch err
                msg = err.message;
                return;
            end
        end
    end
    
end

