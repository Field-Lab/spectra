classdef EditHandler < handle
    %EDITHANDLER object to store and visualize edit commands during cluster edition
    %   The EditHandler class stores internally the sequence of edit options
    %   requested by the user over a dataset in the cluster edition GUI
    %
    %   Edit commands are listed internally as a cell array of strings
    %   And can be viewed in a window with a table listing them
    %
    %   The save button allows to write the edit commands to a file
    %   A prompt will ask user for saving or discarding at the end of edition
    %
    %
    %   Author - Vincent Deo - Stanford University - March 1st 2016
    
    properties(SetAccess = immutable, GetAccess = private)
        dataFolder
        editFilePath
    end
    
    properties(SetAccess = private, GetAccess = public)
        
        allActions
        newActions
        
        filter
        expectedHardClear
        
        editList = cell(0,3) % n x 2 cell array, encoded list of edits requested
        % each 1st col cell is a EditAction enum object
        % each 2nd col cell should be a 1 x k array defining the edit and parameters
        % each 3rd col cell is a 1 x k array giving returned data by backend execution
        displayList = cell(0,1) % n x 1 cell array for display, string info of edit action
        unsavedActions = false
        isUnsaved = false(0,1);
        
        isWindowOn = false
        actionTable = []
        windowHandle = []
        statusBar = []
    end
    
    methods
        % Construtor
        %   Creates the EditHandler object, storage capacities, data table
        %   Input:
        %       datapath: path to the analysis folder, kept here for saving
        function obj = EditHandler(dataFolder)
            obj.dataFolder = dataFolder;
            [~,datasetName,~] = fileparts(dataFolder);
            % Find the manual edit file
            obj.editFilePath = [dataFolder,filesep,datasetName,'.edit.mat'];
            if exist(obj.editFilePath,'file') == 2
                load(obj.editFilePath);
                obj.allActions = manualActions;
                % Find actions not yet consolidated
                lastConsolidate = find(cellfun(@(x) x  == EditAction.CONSOLIDATE,manualActions(:,1),'uni',true),1,'last');
                if numel(lastConsolidate) == 0
                    lastConsolidate = 0;
                end
                obj.newActions = manualActions((lastConsolidate+1):end,:);
            else
                obj.allActions = cell(0,3);
                obj.newActions = cell(0,3);
            end
        end
        
        % function openWindow
        %   opens the GUI window if unexistent or previous one was closed
        function openWindow(obj)
            if obj.isWindowOn
                return
            end
            p = get(groot,'Screensize');
            h = p(4); w = p(3);
            % Pop-up window with title row and type-in text row showing current permutation
            obj.windowHandle = figure( 'Name','Edit Actions Requested', ...
                'MenuBar', 'none', ...
                'Toolbar', 'none', ...
                'NumberTitle', 'off',...
                'Visible', 'on',...
                'OuterPosition', [0.6*w,0.25*h,0.25*w,0.7*h],...
                'deleteFcn',@obj.windowClose);
            supportPanel = uiextras.VBox('Parent',obj.windowHandle,...
                'Position',[0 0 1 1]);
            
            % Cluster data table %
            columnName = {'Action List'};
            columnFormat = {'char'};
            columnEdit = false;
            colWidth =  {321};
            obj.actionTable = uitable(...
                'Parent',supportPanel,...
                'ColumnName', columnName,...
                'ColumnFormat', columnFormat,...
                'ColumnEditable', columnEdit,...
                'Interruptible','off',...
                'Units','norm',...
                'Position',[0 0 1 1],...
                'RowName',[]);
            % Java tweaking to get autoresizing rows for multiline contents
            jscroll = findjobj(obj.actionTable);
            jtable = jscroll.getViewport.getView;
            jtable.setRowAutoResizes(true);
            obj.isWindowOn = true;
            obj.actionTable.ColumnWidth = colWidth;
            obj.updateGUITable();
            
            obj.statusBar = uicontrol(...
                'Parent',supportPanel,...
                'Style','text',...
                'fontsize',10,...
                'String','',...
                'HorizontalAlignment','left');
            supportPanel.Sizes = [-1 24];
        end
        
        % function updateGUITable
        %   updates the uitable contents from the object's listings
        %   if the listing window is open
        function updateGUITable(obj)
            if obj.isWindowOn
                obj.actionTable.Data = obj.displayList; drawnow;
                obj.actionTable.ColumnWidth = cellfun(@plus,obj.actionTable.ColumnWidth,{1},'uni',false); drawnow;
                obj.actionTable.ColumnWidth = cellfun(@plus,obj.actionTable.ColumnWidth,{-1},'uni',false);
            end
        end
        
        % function windowClose()
        %   Callback for when the EditHandler listing window is closed
        %   clears the ui handles and sets the boolean tag to false
        function windowClose(obj,varargin)
            obj.isWindowOn = false;
            obj.actionTable = [];
            obj.windowHandle = [];
        end
            
        % function saveEdits
        %   Writes the edit curently stored in the edit table to a .edit TODO
        %   file in the analysis folder.
        function status = saveEdits(obj)
            if obj.unsavedActions
                if ~obj.isWindowOn
                    obj.openWindow
                else
                    figure(obj.windowHandle);
                end
                choice = questdlg('Do you want to save the pending edits on this electrode?','Save','Save','Clear session','Cancel','Cancel');
                switch choice
                    case 'Save'
                        if obj.expectedHardClear
                            obj.clearHard();
                            obj.expectedHardClear = false;
                        end
                        obj.allActions = [obj.allActions ; obj.editList];
                        obj.newActions = [obj.newActions ; obj.editList];
                        manualActions = obj.allActions;
                        save(obj.editFilePath,'manualActions','-v7.3');
                        obj.unsavedActions = false;
                        obj.isUnsaved(:) = false;
                        status = 0;
                    case 'Clear session'
                        status = 0;
                        obj.clearSession();
                        obj.expectedHardClear = false;
                    case 'Cancel'
                        status = 1;
                end
            else
                status = 0;
            end
        end
        
        % function addAction
        %       add an action and its parameters to the EditHandler storage
        %   Inputs:
        %       action: EditAction object
        %       parameters: cell array describing action parameters
        %       data: backend returned data
        %       mode: unsaved action flag - true unsaved - false saved
        function addAction(obj,action,parameters,data,mode)
            [v,m] = action.checkParameters(parameters);
            if v == 0
                throw(MException('',['EditHandler:addAction - Invalid action parameters: \n',m]));
            end
            obj.editList = [obj.editList ; {action, parameters, data}];
            obj.displayList = [obj.displayList ; obj.genString(action,parameters,data)];
            obj.updateGUITable();
            if mode
                obj.unsavedActions = true;
                obj.isUnsaved = [obj.isUnsaved ; true];
            else
                obj.isUnsaved = [obj.isUnsaved ; false];
            end
        end
        
        % function genString
        %   generates human-readable information string for an action
        %   Inputs:
        %       action: EditAction object
        %       parameters: cell array describing action parameters
        %       data: data cell array resulting from the action execution by the backend
        function str = genString(obj,action,parameters,data)
            validateattributes(action,{'EditAction'},{});
            validateattributes(parameters,{'cell'},{});
            switch action
                case EditAction.KEEP
                    str = ['Elevated status for IDs ', prettyPrint(parameters{1})];
                case EditAction.DELETE
                    str = ['Deletion for IDs ', prettyPrint(parameters{1})];
                case EditAction.MERGE
                    str = ['<html><left />',...
                        'Merge of IDs ',prettyPrint(parameters{1}),'<br />',...
                        'Merge properties:<br />',...
                        '#Spikes: ',num2str(data{1}),' --- Rate: ',num2str(data{2}),...
                        ' Hz --- Contam: ',num2str(data{3}),...
                        '</html>'];
                case EditAction.RECLUSTER
                    str = ['<html><left />',...
                        'Reclustering of IDs ',prettyPrint(parameters{1}),'<br />',...
                        'In ',num2str(parameters{2}),' clusters.<br />',...
                        'With configuration: ''',num2str(parameters{3}),'''<br />',...
                        'Resulting IDs: ',prettyPrint(data{1})];
                case EditAction.SHRINK
                    str = ['<html><left />',...
                        'Shrinking of IDs ',prettyPrint(parameters{1}),'<br />',...
                        'Target fraction: ',num2str(parameters{2}),...
                        ' --- Achieved: ',prettyPrint(data{1}),'<br />',...
                        'Target contamination: ',num2str(parameters{3}),...
                        ' --- Achieved: ',prettyPrint(data{2})];
                otherwise
                    throw(MException('','EditHandler:genString - Unhandled EditAction in switch statement.'));
            end
        end
        
        % function findActionsForElectrode
        %   returns all the new (unconsolidated) actions that apply to the input electrode
        %
        %   Input: el - matlab numbered electrode
        %   Return: actionList - relevant actions
        function actionList = findActionsForElectrode(obj, el)
            IDRange = (15*(el-2)+1):(15*(el-2)+15);
            obj.filter = false(size(obj.newActions,1),1);
            for n = 1:size(obj.newActions,1)
                % All actions have the ID list in params{1}
                if numel(intersect(IDRange,obj.newActions{n,2}{1})) > 0
                    obj.filter(n) = true;
                end
            end
            actionList = obj.newActions(obj.filter,:);
        end
        
        % function removeLastAction
        %   removes the last action in the current actionList
        function removeLastAction(obj)
            obj.editList = obj.editList(1:(end-1),:);
            obj.displayList = obj.displayList(1:(end-1));
            obj.updateGUITable();
        end
        
        % function clearLoad
        %   clears all actions stored by the EditHandler
        %   more of a backend function to flush the window
        function clearLoad(obj)
            obj.editList = cell(0,3);
            obj.displayList = cell(0,1);
            obj.unsavedActions = false;
            obj.isUnsaved = false(0,1);
            obj.updateGUITable();
        end
        
        % function clearSession
        %   clears all edits upon the last saved version
        function clearSession(obj)
            obj.unsavedActions = false;
            obj.expectedHardClear = false;
            obj.editList = obj.newActions(obj.filter,:);
            obj.displayList = cellfun(@(t,u,v) obj.genString(t,u,v),obj.editList(:,1),obj.editList(:,2),obj.editList(:,3),'uni',false);
            obj.isUnsaved = false(size(obj.editList,1),1);
            obj.updateGUITable();
        end
        
        % function requestClearHard()
        %   request a hard clearing that will be applied only at the next save prompt.
        function requestClearHard(obj)
            obj.expectedHardClear = true;
            obj.unsavedActions = true;
        end
        
        % function closeWindow
        %   close the editHandlerWindow
        function closeWindow(obj)
            if obj.isWindowOn
                delete(obj.windowHandle);
                obj.isWindowOn = false;
            end
        end
    end
    
    methods(Access = private)
        % function clearHard
        %   clear actions for the electrode and removes them from the
        %   allActions/actionNew storage
        %   and for that reason they will be deleted from the .edit.mat file
        %   at next save, even if they were there before
        %
        %   Is private and should only be called when actually saving.
        %   because if the user choses to discard session at the sve prompt,
        %   then the hard delete should not be performed at all
        function clearHard(obj)
            obj.newActions(obj.filter,:) = [];
            obj.allActions([false(size(obj.allActions,1)-numel(obj.filter),1),obj.filter],:) = [];
            obj.filter(obj.filter) = [];
        end
            
    end % PRIVATE methods
    
end

