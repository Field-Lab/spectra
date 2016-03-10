classdef EditAction
    %EDITACTION enumeration class describing possible edition actions
    %   This class specifies the possible edition action and
    %   Their required parameters are defined by the method checkParameters
    %
    %   See ClusterEditBackend for how actions are applied
    %   and ClusterEditGUI for how they are called
    
    enumeration
        DEBUG
        DEBUG_2
        
        NO_REMOVE
        MERGE
        NO_MERGE
        SHRINK
        RECLUSTER
    end
    
    methods
        % function checkParameters
        %   Checks if the parameter cell array params has a valid structure for action obj
        %   This function is the definition of what are the required parameters for an action and
        %   their format.
        %
        %   Input:
        %       params: parameter cell array
        %   Returns:
        %       valid: boolean = params is a valid parameters array for action obj
        %       msg: error message in case params is invalid.
        %
        % --- Detailed expected parameters ---
        %   DEBUG:      {nonempty row vector}
        %   DEBUG_2:    {nonempty row vector , string}
        %   NO_REMOVE:     {nonempty row vector: neuron IDs}
        %
        function [valid,msg] = checkParameters(obj,params)
            valid = false;
            msg = '';
            unhandled = false;
            try
                switch obj
                    case EditAction.DEBUG
                        validateattributes(params,{'cell'},{'size',[1 1]},'','parameter cell array');
                        validateattributes(params{1},{'numeric'},{'row','nonempty'},'','DEBUG random nonempty list',1);
                    case EditAction.DEBUG_2
                        validateattributes(params,{'cell'},{'size',[2 1]},'','parameter cell array');
                        validateattributes(params{1},{'numeric'},{'row','nonempty'},'','DEBUG_2 random nonempty list',1);
                        validateattributes(params{2},{'char'},{'row','nonempty'},'','DEBUG_2 nonempty string',2);
                    case EditAction.NO_REMOVE
                        validateattributes(params,{'cell'},{'size',[1 1]},'','parameter cell array');
                        validateattributes(params{1},{'numeric'},{'row','nonempty'},'','List of IDs to elevate',1);
                    otherwise
                        unhandled = true; % report throw out of try block
                end
                valid = true;
            catch err
                msg = err.message;
                return;
            end
            if unhandled
                throw(MException('','EditAction:checkParameters - Unhandled EditAction in switch statement.'));
            end
        end
        
        function s = getTooltipString(obj)
            switch obj
                case EditAction.DEBUG
                    s = 'Debug action, does nothing.';
                case EditAction.DEBUG_2
                    s = 'Debug action, does nothing.';
                case EditAction.NO_REMOVE
                    s = 'Elevate the neurons statuses. They won''t be affected by duplicate removal until another calculation is ran.';
                case EditAction.MERGE
                    s = 'Merge together the cluster selection.';
                case EditAction.NO_MERGE
                    s = 'Forbid selected clusters to be merged with any other cluster.';
                case EditAction.SHRINK
                    s = ['Shrink a cluster to reach a given percentage of initial size,',...
                        'following a topologically homogeneous process that conserves shape.'];
                case EditAction.RECLUSTER
                    s = 'Perform spectral clustering form scratch on the selected clusters, with fixed or automatic number of clusters.';
                otherwise
                    throw(MException('','EditAction:getTooltipString - Unhandled EditAction in switch statement.'));
            end
        end
    end
    
end

