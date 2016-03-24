classdef EditAction
    %EDITACTION enumeration class describing possible edition actions
    %   This class specifies the possible edition action and
    %   Their required parameters are defined by the method checkParameters
    %
    %   See ClusterEditBackend for how actions are applied
    %   and ClusterEditGUI for how they are called
    properties
        isManual % Tag filtering between manual actions and automatic action.
        % Manual actions only should generate a button and be user accessible in the GUI
    end
    
    enumeration
        CONSOLIDATE  (0) % CONSOLIDATE action is a marker for incrementally edited neurons files.
        % It should only be added automatically when neurons files are rebuilt.
        
        AUTO_RM      (0)
        AUTO_MERGE   (0)
        AUTO_RM_DUP  (0)
        
        ELEVATE      (1)
        MERGE        (1)
        SHRINK       (1)
        RECLUSTER    (1)
    end
     
    methods
        function c = EditAction(x)
            c.isManual = x;
        end
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
        %   ELEVATE:      {nonempty row vector: neuron IDs}
        %   MERGE:          {nonempty, nonsingleton row vector: neuron IDs}
        %   RECLUSTER:      {nonempty row vector: neuron IDs,
        %                       positive integer scalar: number of clusters to make. 0 for auto,
        %                       string: configuration tag to use for spectral clustering }
        %   SHRINK:         {nonempty row vector: neuron IDs,
        %                       scalar in [0,1]: fraction of spikes to keep at most,
        %                       scalar >= 0: Contamination to achieve at most}
        function [valid,msg] = checkParameters(obj,params)
            valid = false;
            msg = '';
            unhandled = false;
            try
                switch obj
                    case EditAction.ELEVATE
                        validateattributes(params,{'cell'},{'size',[1 1]},'','parameter cell array');
                        validateattributes(params{1},{'numeric'},{'row','nonempty'},'','List of IDs to elevate',1);
                    case EditAction.MERGE
                        validateattributes(params,{'cell'},{'size',[1 1]},'','parameter cell array');
                        validateattributes(params{1},{'numeric'},{'row','nonempty'},'','List of IDs to merge',1);
                        validateattributes(numel(params{1}),{'numeric'},{'scalar','>=',2},'','List of IDs to merge (numel)',1);
                    case EditAction.RECLUSTER
                        validateattributes(params,{'cell'},{'size',[3 1]},'','parameter cell array');
                        validateattributes(params{1},{'numeric'},{'row','nonempty'},'','List of IDs to recluster',1);
                        validateattributes(params{2},{'numeric'},{'scalar','integer','>=',0},'','Number of clusters to apply',2);
                        validateattributes(params{3},{'char'},{},'','Configuration tag',3);
                    case EditAction.SHRINK
                        validateattributes(params,{'cell'},{'size',[3 1]},'','parameter cell array');
                        validateattributes(params{1},{'numeric'},{'row','nonempty'},'','List of IDs to recluster',1);
                        validateattributes(params{2},{'numeric'},{'scalar','>=',0,'<=',1},'','Target fraction of original spikes',2);
                        validateattributes(params{3},{'numeric'},{'scalar','>=',0},'','Target Contamination',3);
                        if params{3} == Inf && params{2} == 1
                            throw(MException('','Will not shrink with values (1,Inf). Equivalent to nothing.'));
                        end
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
                case EditAction.ELEVATE
                    s = 'Elevate the neurons statuses. They won''t be affected by duplicate removal or merges until another calculation is ran.';
                case EditAction.MERGE
                    s = 'Merge together the cluster selection.';
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
    
    methods(Static)
        % function niceDisplay
        %   prettyprints an N * 3 array of editAction,parameters,data to the console
        %   Common idea with editHandler's genString.
        %   Input:
        %       actionList : N * 3 cell array
        %           [{ EditAction , { parameter }, {data} } ; ... ]
        function s = toStringList(actionList)
            nActions = size(actionList,1);
            s = cell(nActions,1);
            for n = 1:nActions
                action = actionList{n,1};
                params = actionList{n,2};
                data = actionList{n,3};
                
                
                a = sprintf('%s',action.char);
                if numel(params) > 0
                    pString = cellfun(@prettyPrint,params,'un',false);
                else
                    pString = {''};
                end
                p = sprintf('%s    ',pString{:});
                d = evalc('data');
                s{n} = sprintf('%s\n%s\n%s',a,p,d(15:(end-2)));
                % Just putting this here to encourage a future action specific display
                %{
                switch action
                    case EditAction.CONSOLIDATE
                    case EditAction.AUTO_RM
                    case EditAction.AUTO_MERGE
                    case EditAction.AUTO_RM_DUP
                    case EditAction.ELEVATE
                    case EditAction.MERGE
                    case EditAction.SHRINK
                    case EditAction.RECLUSTER
                    otherwise
                        throw(MException('','EditAction:toStringList - Unhandled EditAction in switch statement.'));
                end
                %}
            end
            
            s = sprintf('%s\n--------------------\n',s{:});
        end
    end % Static Methods
    
end

