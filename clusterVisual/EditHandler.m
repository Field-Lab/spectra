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
    %   To apply edition changes, the user should make a call to the
    %   matlab TODO function, which will apply the provided changes,
    %   then repeat all the non-interfering duplicate removal steps listed in
    %   the .clean.mat file
    %   See the documentation of TODO function for how changes will be applied
    %   to the dataset
    %
    %   Author - Vincent Deo - Stanford University - March 1st 2016
    
    properties(SetAccess = immutable, GetAccess = private)
        dataFolder
    end
    
    properties(SetAccess = private, GetAccess = private)
        editList = cell(0,2) % n x 2 cell array, encoded list of edits requested
        % each 1st col cell is a EditAction enum object
        % each 2nd col cell should be a 1 x k array defining the edit and parameters
        displayList = cell(0,1) % n x 1 cell array for display, string info of 
        
        isWindowOn = false
        actionTable = []
        windowHandle = []
        
        lastDel % deletion recovery
    end
    
    methods
        % Construtor
        %   Creates the EditHandler object, storage capacities, data table
        %   Input:
        %       datapath: path to the analysis folder, kept here for saving
        function obj = EditHandler(dataFolder)
            obj.dataFolder = dataFolder;
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
                'OuterPosition', [0.7*w,0.25*h,0.25*w,0.7*h],...
                'deleteFcn',@obj.windowClose);
            supportPanel = uipanel('Parent',obj.windowHandle,...
                'Position',[0 0 1 1]);
            
            % Cluster data table %
            columnName = {'Action List','Del'};
            columnFormat = {'char','logical'};
            columnEdit = [false, true];
            colWidth =  {0.25*w-47,26};
            obj.actionTable = uitable(...
                'Parent',supportPanel,...
                'ColumnName', columnName,...
                'ColumnFormat', columnFormat,...
                'ColumnEditable', columnEdit,...
                'Units','norm',...
                'Position',[0 0 1 1],...
                'RowName',[],...
                'CellEditCallback',@obj.delRow);
            % Java tweaking to get sortable columns
%             jscrollpane = findjobj(obj.actionTable);
%             jtable = jscrollpane.getViewport.getView;
%             jtable.setSortable(true);
%             jtable.setAutoResort(true);
%             jtable.setMultiColumnSortable(true);
%             jtable.setPreserveSelectionsAfterSorting(true);
            obj.actionTable.ColumnWidth = colWidth;
            % End cluster data table %
            obj.isWindowOn = true;
            obj.updateGUITable();
        end
        
        % function delRow
        %   removes an action from the object's action listing.
        %   Is the callback of the uitable (second column); when the checkbox is clicked,
        %   the row's action is deleted
        function delRow(obj,source,callbackdata)
            obj.lastDel.pos = callbackdata.Indices(1);
            obj.lastDel.display = obj.displayList(callbackdata.Indices(1),:);
            obj.lastDel.edit = obj.editList(callbackdata.Indices(1),:);
            obj.displayList(callbackdata.Indices(1),:) = [];
            obj.editList(callbackdata.Indices(1),:) = [];
            
            obj.updateGUITable;
        end
        
        % function restoreLastDel
        %   restores the latest deleted action back into the action listings
        %   Does nothing if there is no action to restore
        %   Cannot be called twice in a row
        function restoreLastDel(obj)
            if numel(obj.lastDel) > 0
                obj.editList = [obj.editList(1:(obj.lastDel.pos-1),:);...
                    obj.lastDel.edit;...
                    obj.editList((obj.lastDel.pos):end,:)];
                obj.displayList = [obj.displayList(1:(obj.lastDel.pos-1),:);...
                    obj.lastDel.display;...
                    obj.displayList((obj.lastDel.pos):end,:)];
                
                obj.lastDel = [];
                obj.updateGUITable;
            end
        end
        
        % function updateGUITable
        %   updates the uitable contents from the object's listings
        %   if the listing window is open
        function updateGUITable(obj)
            if obj.isWindowOn
                obj.actionTable.Data = [obj.displayList,num2cell(false(numel(obj.displayList),1))];
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
        function saveEdits(obj)
            if ~obj.isWindowOn
                obj.openWindow
            else
                figure(obj.windowHandle);
            end
            path = [obj.dataFolder,filesep,'edits.mat'];
            if exist(path)
                choice = questdlg('An edit file was found...','Save edit commands','Overwrite','Append','Cancel','Cancel');
            else
                choice = questdlg('Save edits ?','Save edit commands','Save','Cancel','Save');
            end
            editList = obj.editList;
            switch choice
                case 'Save'
                    save(path,'editList','-v7.3');
                case 'Overwrite'
                    save(path,'editList','-v7.3');
                case 'Append'
                    save(path,'editList','-append');
            end
        end
        
        % function addAction
        %       add an action and its parameters to the EditHandler storage
        %   Inputs:
        %       action: EditAction object
        %       parameters: cell array describing action parameters
        function addAction(obj,action,parameters)
            validateattributes(action,{'EditAction'},{});
            validateattributes(parameters,{'cell'},{});
            obj.editList = [obj.editList ; {action, parameters}];
            obj.displayList = [obj.displayList ; obj.genString(action,parameters)];
            obj.updateGUITable();
        end
        
        % function genString
        %   generates human-readable information string for an action
        %   Inputs:
        %       action: EditAction object
        %       parameters: cell array describing action parameters
        function str = genString(obj,action,parameters)
            validateattributes(action,{'EditAction'},{});
            validateattributes(parameters,{'cell'},{});
            switch action
                case EditAction.DEBUG
                    str = sprintf('DEBUG Action with %u parameters',numel(parameters));
                case EditAction.DEBUG_2
                    str = sprintf('DEBUG_2 Action with %u parameters',numel(parameters));
            end
        end
        
        % function clearActions
        %   clears all actions stored by the EditHandler
        function clearActions(obj)
            obj.editList = cell(0,2);
            obj.displayList = cell(0,1);
            obj.updateGUITable();
        end
    end
    
end

