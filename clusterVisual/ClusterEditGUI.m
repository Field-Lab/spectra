function varargout = ClusterEditGUI(datasetFolder,varargin)
    %CLUSTEREDITGUI
    % Manages all components of the GUI for the cluster visualizer
    % Placement of all graphical boxes
    % Callbacks for all user interaction
    % Which makes quite a large file
    %
    % Input argument:
    %       datasetFolder: path to analysis folder
    %           Must mandatorily contain a neurons.mat and a prj.mat
    %           Optionally, .clean.mat, .ei, .sta, and .classification.txt
    %           Will be read and provide additional information
    %       OPTIONAL - backEndHandle
    %           Provides a shortcut to restart figure without building
    %           a new data handling object of class ClusterEditBackEnd
    %
    % Outputs:
    %       OPTIONAL - frontEndHandle (1st)
    %           Handle to main GUI figure
    %       OPTIONAL - backEndHandle (2nd)
    %           handle to clusterEditBackend data handler
    %
    % Author - Vincent Deo - Stanford University - Feb 29th 2016
    
    %% VERSION NUMBER %%
    % Increment at each master merge
    version = '0.2';
    
    %% Main frontend and backend startup %%
    narginchk(1,2); % Existing back-end may be passed in varargin{1}
    if nargin == 1 % New backend
        backEndHandle = ClusterEditBackend(datasetFolder);
    else % Existing backend - overrides first argument datasetFolder
        backEndHandle = varargin{1};
        datasetFolder = backEndHandle.analysisPath;
    end
    
    % Attach edit handler
    editHandler = EditHandler(datasetFolder);
    
    % Front-end main figure
    frontEndHandle = figure( 'Name',sprintf('Cluster Editor %s',version), ...
        'MenuBar', 'none', ...
        'Toolbar', 'figure', ...
        'NumberTitle', 'off',...
        'Visible', 'off',...
        'OuterPosition',get(groot,'Screensize') + [60 60 -120 -90],...
        'DeleteFcn',@closeWinCallback);
    % Remove buttons - we only keep figure movement handling
    [~,~] = customizeTools(frontEndHandle);
       
    %% Configuration %%
    spacerWidth = 10; % generic spacer size (pixel; but hardcoded values are often used)
    displayPoints = 8000; % plot display
    clusterColors = zeros(0,3); % Empty global initializer for display colors
    
    %% GUI layout %%
    % Main box
    mainLayout = uiextras.HBoxFlex(...
        'Parent',frontEndHandle,...
        'Position',[0 0 1 1],...
        'Spacing',2);
    
    % Left VBox - Menu, title, table, EI/STA displays
    leftColumns = uiextras.VBox(...
        'Parent',mainLayout,...
        'Spacing',spacerWidth);
    
    % Right VBox - PCs, ACF, Spike rate plots
    graphLayout = uiextras.VBoxFlex(...
        'Parent',mainLayout,...
        'Spacing',spacerWidth);
    
    mainLayout.Sizes = [-1 -1]; % Finish mainLayout
    
    % 3D plot box - HBox will allow to add a vertical button strip to the side
    % of the 3D plot
    graph3DBox = uiextras.HBox(...
        'Parent',graphLayout,...
        'Spacing',3);
    
    % Panel for 3D box axes (always put a panel to encapsulate axes)
    supportPanel3D = uipanel('Parent',graph3DBox);
    
    % PC ordering
    PC = 1:5;
    
    % 3D plot initialization %
    PC123Plots = {}; % handle to globalize 3D scatter plot access
    plot3D = axes(...
        'Parent',supportPanel3D,...
        'ClippingStyle','rectangle',...
        'View',[45,15],...
        'XGrid','on','YGrid','on','ZGrid','on',...
        'DataAspectRatio',[1 1 1],...
        'XLimMode','manual','YLimMode','manual','ZLimMode','manual');
    plot3D.Title.String = sprintf('Principal Components %u-%u-%u',PC(1),PC(2),PC(3));
    plot3D.XLabel.String = sprintf('PC %u',PC(1));
    plot3D.YLabel.String = sprintf('PC %u',PC(2));
    plot3D.ZLabel.String = sprintf('PC %u',PC(3));
    axis(plot3D,'tight');
    % End 3D plot initialization
    
    % Rightmost column of edition buttons
    editCol = uiextras.VBox(...
        'Parent',graph3DBox,...
        'Spacing',2);
    uiextras.Empty('Parent',graph3DBox); % empty colon at right of toolstrip, to recover spacing at edge of screen
    editEnum = enumeration('EditAction');
    openEditHandler = uicontrol('Parent',editCol,'Style','pushbutton',...
        'String','-SHOW LIST-',...
        'TooltipString','Show the edit viewer edit list',...
        'Interruptible','off','BusyAction','cancel',...
        'fontsize',9,'callback',@openEditHandlerCallback);
    saveEditActions = uicontrol('Parent',editCol,'Style','pushbutton',...
        'String','-SAVE EDITS-',...
        'TooltipString','Save the edit actions to the .edit.mat file',...
        'Interruptible','off','BusyAction','cancel',...
        'fontsize',9,'callback',@saveEditsCallback);
    clearEditActions = uicontrol('Parent',editCol,'Style','pushbutton',...
        'String','-CLEAR-',...
        'TooltipString','Clear all edit actions from this session on this electrode.',...
        'Interruptible','off','BusyAction','cancel',...
        'fontsize',9,'callback',@clearEditsCallback);
    deleteEditActions = uicontrol('Parent',editCol,'Style','pushbutton',...
        'String','-DEL ALL-',...
        'TooltipString','Clear all edit actions on this electrode - including those stored in .edit.mat.',...
        'Interruptible','off','BusyAction','cancel',...
        'fontsize',9,'callback',@clearEditsCallback);
    undoLastEditButton = uicontrol('Parent',editCol,'Style','pushbutton',...
        'String','-UNDO-',...
        'TooltipString','Undo latest action (from this session only) on this electrode.',...
        'Interruptible','off','BusyAction','cancel',...
        'fontsize',9,'callback',@undoLastEdit);
    uiextras.Empty('Parent',editCol,'background','g');
    actionButtons = cell(numel(editEnum),1);
    for actionNum = 1:numel(editEnum)
        action = editEnum(actionNum);
        if action.isManual
            str = action.char;
            str = [str(1),lower(str(2:end))];
            actionButtons{actionNum} = uicontrol(...
                'Parent',editCol,...
                'Style', 'pushbutton',...
                'String',str,...
                'TooltipString',action.getTooltipString,...
                'Interruptible','off','BusyAction','cancel',...
                'fontsize',9,...
                'callback',{@editCallback,action});
        end
    end
    uiextras.Empty('Parent',editCol,'background','g');
    editCol.Sizes = [24 * (ones(1,6)), repmat(24,1,nnz(arrayfun(@(x) x.isManual,editEnum))),-1];
    
    graph3DBox.Sizes = [-1 75 0]; % Finish graph3DBox
    
    % HBox wrapper for [[Sp rate; ACF], PC45]
    bottomGraphs = uiextras.HBoxFlex(...
        'Parent',graphLayout,...
        'Spacing',spacerWidth);
    graphLayout.Sizes = [-3 -2]; % Finish graphLayout Vbox
    % VBox wrapper for [Sp rate ; ACF]
    bottomGraphsLeftBox = uiextras.VBox(... 
        'Parent',bottomGraphs,...
        'Spacing',0);
    % Panel wrapper for PC45
    bottomGraphsRightPanel = uipanel(...
        'Parent',bottomGraphs);
    bottomGraphs.Sizes = [-1 -1]; % Finish bottom right graph panel
    
    % Spike rate panel
    ratePanel = uipanel(...
        'Parent',bottomGraphsLeftBox);
    % Initialize Spike Rate plot %
    ratePlots = {}; % global handle
    ratePlot = axes(...
        'Parent',ratePanel,...
        'XGrid','on','YGrid','on',...
        'OuterPosition',[0 0 1 1]);
    ratePlot.Title.String = 'Spike rates'; ratePlot.Title.FontSize = 9;
    ratePlot.XLabel.String = 'Time (sec)'; ratePlot.XLabel.FontSize = 9;
    ratePlot.YLabel.String = 'Spike rate (Hz)'; ratePlot.YLabel.FontSize = 9;
    ratePlot.XLim = [0, backEndHandle.nSamples / 20000 + 1];
    % end spike rate plot %
    
    % Panel to hold ACF plot
    ACFPanel = uipanel(...
        'Parent',bottomGraphsLeftBox);
    % Initialize ACF plot %
    ACFPlots = {}; % global handle
    ACFPlot = axes(...
        'Parent',ACFPanel,...
        'XGrid','on','YGrid','on',...
        'OuterPosition',[0 0 1 1]);
    ACFPlot.Title.String = 'ACF'; ACFPlot.Title.FontSize = 9;
    ACFPlot.XLabel.String = 'Time \Delta (msec)'; ACFPlot.XLabel.FontSize = 9;
    ACFPlot.YLabel.String = 'Autocorr. (pair fraction/msec \Delta)'; ACFPlot.YLabel.FontSize = 9;
    ACFPlot.XLim = [0, 101];
    % end ACF plot %
    
    bottomGraphsLeftBox.Sizes = [-1 -1]; % Finish [Sp rate ; ACF] VBox
    
    % Initialize PC45 plot %
    PC45Plots = {}; % global handle
    PC45Plot = axes(...
        'Parent',bottomGraphsRightPanel,...
        'XGrid','on','YGrid','on',...
        'OuterPosition',[0 0 1 1],...
        'XLimMode','manual','YLimMode','manual','ZLimMode','manual');
    PC45Plot.Title.String = sprintf('Principal Components %u-%u',PC(4),PC(5));
    PC45Plot.XLabel.String = sprintf('PC %u',PC(4));
    PC45Plot.YLabel.String = sprintf('PC %u',PC(5));
    axis(PC45Plot,'tight');
    % End PC45 Plot %
    
    % Contents of the left half VBox
    % Dataset name - title
    datasetName = uicontrol(...
        'Parent',leftColumns,...
        'Style','text',...
        'fontsize',12,...
        'String',datasetFolder);
    % VBox encapsulating control features (button rows and table)
    menu = uiextras.VBox(...
        'Parent',leftColumns,...
        'Spacing',0);
    % Bottom left box - wrapping EI, STA, EI dist matrix displays
    eistaBox = uiextras.HBox('Parent',leftColumns,...
        'Spacing',0);
    % Status bar at bottom of left column
    statusBar = uicontrol(...
        'Parent',leftColumns,...
        'Style','text',...
        'fontsize',10,...
        'String','Welcome to Cluster Editor 0.0 !',...
        'HorizontalAlignment','left');
    
    % Pass status bar handle to backend
    backEndHandle.statusBarHandle = statusBar;
    
    leftColumns.Sizes = [24 -1 0 24]; % Finish general layout of left half
    
    % Menu VBox %
    % First row - load buttons, text boxes
    loadRow = uiextras.HBox(...
        'Parent',menu,...
        'Spacing',6,...
        'Padding',2);
    % Text box, electrode number to load
    elNumberBox = uicontrol(...
        'Parent',loadRow,...
        'Style', 'edit',...
        'Max',1,'Min',1,...
        'String','El#',...
        'TooltipString','Type in an electrode number and hit return to load.',...
        'fontsize',10,...
        'callback',@loadButtonCallback);
    % increment/decrement buttons box
    ppmmButtonBox = uiextras.VBox(...
        'Parent',loadRow,...
        'Spacing',0,...
        'Padding',0);
    % Electrode increment button
    ppButton = uicontrol(...
        'Parent',ppmmButtonBox,...
        'Style', 'pushbutton',...
        'String',char(hex2dec('25B2')),...
        'TooltipString','Increment electrode # and load.',...
        'fontsize',6,...
        'callback',@loadButtonCallback);
    % Electrode decrement button
    mmButton = uicontrol(...
        'Parent',ppmmButtonBox,...
        'Style', 'pushbutton',...
        'String',char(hex2dec('25BC')),...
        'TooltipString','Decrement electrode # and load.',...
        'fontsize',6,...
        'callback',@loadButtonCallback);
    ppmmButtonBox.Sizes = [-1 -1]; % Finish increment/decrement button box
    % Button actively loading typed electrode
    loadButton = uicontrol(...
        'Parent',loadRow,...
        'Style', 'pushbutton',...
        'String','Load',...
        'TooltipString','Click to load the electrode # typed in the box.',...
        'fontsize',10,...
        'callback',@loadButtonCallback);
    % Text box, cluster ID to seek and load
    clustNumberBox = uicontrol(...
        'Parent',loadRow,...
        'Style', 'edit',...
        'Max',1,'Min',1,...
        'String','ID#',...
        'TooltipString','Type in a neuron ID and hit return to load its electrode.',...
        'fontsize',10,...
        'callback',@loadClustButtonCallback);
    % Button actively loading ID typed in text box
    loadClusterButton = uicontrol(...
        'Parent',loadRow,...
        'Style', 'pushbutton',...
        'String','Load',...
        'TooltipString','Click to load neuron ID typed in the box',...
        'fontsize',10,...
        'callback',@loadClustButtonCallback);
    % Filler
    fillerLoadRow = uiextras.Empty(...
        'Parent',loadRow);
    loadRow.Sizes = [ 34, 20, 50, 45, 60, -1]; % finish loadRow
    
    % selector row w/ "(un)select all", re-scatter and rescale axes buttons
    selectorRow = uiextras.HBox(...
        'Parent',menu,...
        'Spacing',6,...
        'Padding',2);
    % Select all clusters button
    allButton = uicontrol(...
        'Parent',selectorRow,...
        'Style', 'pushbutton',...
        'String','Select All',...
        'TooltipString','Display of all the neurons in the table below.',...
        'fontsize',10,...
        'callback',@selectorCallback);
    % Unselect all clusters button
    noneButton = uicontrol(...
        'Parent',selectorRow,...
        'Style', 'pushbutton',...
        'String','Unselect All',...
        'TooltipString','Hide all the neurons',...
        'fontsize',10,...
        'callback',@selectorCallback);
    % filler
    fillerSelectorRow = uiextras.Empty(...
        'Parent',selectorRow);
    % Flushed right : resample displayed points button
    reScatterButton = uicontrol(...
        'Parent',selectorRow,...
        'Style', 'pushbutton',...
        'String','Re-scatter',...
        'TooltipString','Redo random sampling of PC points displayed.',...
        'fontsize',10,...
        'callback',@refreshGraphics);
    % Automatic PC plot axis rescaling
    autoscaleButton = uicontrol(...
        'Parent',selectorRow,...
        'Style', 'pushbutton',...
        'String','Rescale Axes',...
        'TooltipString','Autoscale axes of 3D and PC 4-5 plots.',...
        'fontsize',10,...
        'callback',@autoScalePCPlots);
    switchPCOrderButton = uicontrol(...
        'Parent',selectorRow,...
        'Style', 'pushbutton',...
        'String','Change PCs',...
        'TooltipString','Change PC on axes of 3D and 2D PC plots.',...
        'fontsize',10,...
        'callback',@switchPC);
    selectorRow.Sizes = [75 90 -1 75 95 80]; % finish selectorRow
    
    % Information row on selected electrode
    infoRow = uicontrol(...
        'Parent',menu,...
        'Style','text',...
        'fontsize',10,...
        'String','',...
        'HorizontalAlignment','center');
    
    % Cluster data table %
    columnName = {'ID','','Disp','Status','#Spikes','Rate (Hz)','Contam','Classification'};
    columnFormat = {'numeric','char','logical','char','numeric','numeric','numeric','char'};
    columnEdit = [false, false, true, false, false, false, false, false];
    px = 10; colWidth =  {px*4, px*4, px*3, px*12, px*5, px*6, px*6, px*10};
    d = cell(0,8);
    clustMgmt = uitable(...
        'Parent',menu,...
        'ColumnName', columnName,...
        'ColumnFormat', columnFormat,...
        'ColumnEditable', columnEdit,...
        'RowName',[],...
        'CellEditCallback',@tableEditCallback);
    % Java tweaking to get sortable columns
    jscrollpane = findjobj(clustMgmt);
    jtable = jscrollpane.getViewport.getView;
    jtable.setSortable(true);
    jtable.setAutoResort(true);
    jtable.setColumnAutoResizable(true);
    jtable.setMultiColumnSortable(true);
    jtable.setPreserveSelectionsAfterSorting(true);
    clustMgmt.ColumnWidth = colWidth;
    clustMgmt.Data = d;
    % End cluster data table %
    
    menu.Sizes = [24 24 24 -1]; % Finish left menu
    
    % Initialize EI and STA displays %
    % Left part of bottom left panel
    eiBox = uiextras.VBox(...
        'Parent',eistaBox,...
        'Padding',0,'Spacing',0);
    % EI info at top of eiBox
    eiInfoString = uicontrol(...
        'Parent',eiBox,...
        'Style','text',...
        'fontsize',12,...
        'String','EI Display',...
        'Visible','on');
    eiInfoString.Units = 'Norm';
    % Panel for ei display
    eiPanel = uipanel('Parent',eiBox,...
        'BorderType','none');
    % Option 1: There is an EI - Java panel wrapper to use vision JPanel
    eiJPanel = []; eiHPanel = []; % Globalize the java panels
    % Option 2: There is no EI. Axes to display a graph instead. Unused.
    bottomLeftAxes = axes('Parent',eiPanel,...
        'Visible','off');
    eiBox.Sizes = [26 -1]; % Finish EI display box
    
    % Right part of bottom left panel
    % Option 1: STA to display: Panel and JPanel global variable support
    staPanel = uipanel('Parent',eistaBox,...
        'BorderType','none');
    staJPanel = []; staHPanel = []; % Globalize the java panels
    % String in place of STA for single selected neuron but no STA
    noSTAString = uicontrol(...
        'Parent',staPanel,...
        'Style','text',...
        'fontsize',12,...
        'String','No STA for this neuron',...
        'Visible','off');
    noSTAString.Units = 'Norm';
    noSTAString.Position = [0 0 1 1];
    
    % Option 2: no STA to display - show EI distance matrix
    % VBox taking the place of the STA panel
    bottomRightAxesBox = uiextras.VBox(...
        'Parent',staPanel,'Visible','off');
    % Info/button row on top of the EI matrix display
    bottomRightControlRow = uiextras.HBox(...
        'Parent',bottomRightAxesBox,...
        'Spacing',6);
    % Title string in info row
    EIMatTitleString = uicontrol(...
        'Parent',bottomRightControlRow,...
        'Style','text',...
        'fontsize',12,...
        'String','Pairwise EI Distance   ',...
        'Visible','on');
    % Switch to default permutation button
    defaultPerm = uicontrol(...
        'Parent',bottomRightControlRow,...
        'String','Dflt',...
        'callback',@defaultPermutationCallback);
    % Switch to optimal permutation button
    optPerm = uicontrol(...
        'Parent',bottomRightControlRow,...
        'String','Optml',...
        'callback',@optimalPermutationCallback);
    % Switch to custom permutation button
    customPerm = uicontrol(...
        'Parent',bottomRightControlRow,...
        'String','Custom ...',...
        'callback',@customPermutationCallback);
    currentPermutation = zeros(1,0); % Globalize the permutation of clusters to display
    bottomRightControlRow.Sizes = [-1 45 45 55]; % finish controller row of EI mat display
    % Support panel for EI distance mat axes
    bottomRightAxesPanel = uipanel(...
        'Parent',bottomRightAxesBox,...
        'Visible','on');
    bottomRightAxesBox.Sizes = [26 -1]; % Finish ei mat control row + axes panel box
    % Axes of EI distance matrix
    bottomRightAxes = axes('Parent',bottomRightAxesPanel,...
        'Visible','on');
    
    eistaBox.Sizes = [-1 -1]; % Finish the bottom left part
    
    
    
    %% ---------------------------------------------------------- %%
    % CALLBACKS AND SUBFUNCTIONS
    %---------------------------------------------------------------
    
    % Lock for electrode loader
    % will prevent loader to be interrupted (matlab is actually pretty mono
    % threaded and allows this) to not mess backend and to not uselessly
    % load/display data if another electrode is requested before finishing
    myGlobalLoadInterrupt = 0;
    
    % function loadButtonCallback
    %   Callback of electrode loader interfaces:
    %   electrode textBox, electrode load button,
    %   increment/decrement electrode button
    %
    %   Updates textboxes contents
    %   Requests load to backend
    %   Updates data and graphics
    %
    %   Gracefully handles multiples fast request (multiple clicks on arrow buttons) 
    function loadButtonCallback(source,callbackdata)
        myGlobalLoadInterrupt = 1; % Grab "lock"
        
        e = str2num(elNumberBox.String) + 1; % String is in JAVA numbering, e is in matlab numbering.
        if numel(e) ~= 1 || isnan(e) % Check invalid format
            statusBar.String = sprintf('Requested electrode %s not a number.',elNumberBox.String);
            return;
        end
        if source == ppButton % Incrementer button
            elNumberBox.String = num2str(e);
            e = e+1;
        end
        if source == mmButton % Decrementer button
            elNumberBox.String = num2str(e-2);
            e = e-1;
        end
        statusBar.String = sprintf('Loading electrode %s...',elNumberBox.String);
        
        % Grab the interruption queue if there is one.
        pause(.5); % Give us time to get interrupted (avoid redudancy queuing in case of frantic clicking on ++/-- buttons)
        % If interrupted, do not resume (opposite of matlab default behavior)
        if ~myGlobalLoadInterrupt
            return;
        end
        source.Interruptible = 'off'; % Disable interruption
        
        % TODO: EDIT HANDLER - check for required saving
        status = editHandler.saveEdits();
        if status ~= 0 % Failed save (basically user clicked Cancel)
            return;
        end
        editHandler.clearLoad();
        
        status = backEndHandle.loadEl(e); % Load in backend
        if status ~= 0 % Failed load (not a valid electrode number or already loaded)
            return;
        end
        
        actions = editHandler.findActionsForElectrode(e);
        for i = 1:size(actions,1)
            [~,data] = backEndHandle.localApplyAction(actions{i,1},actions{i,2},actions{i,3});
            editHandler.addAction(actions{i,1},actions{i,2},data,false);
        end
        if numel(actions) > 0
            editHandler.openWindow();
        end
        
        % Allow for interruption again - before refreshing data/graphics
        % Grab interruption queue
        source.Interruptible = 'on';
        drawnow;
        if ~myGlobalLoadInterrupt % If interrupted, do not resume
            return;
        end
        
        % Refresh data and graphics
        statusBar.String = 'Refreshing View...'; drawnow;
        refreshView();
        statusBar.String = sprintf('Displaying electrode %u.',e-1);
        % Reset string displayed in cluster text box if this was a load by electrode
        if source == loadButton || source == elNumberBox || source == ppButton || source == mmButton
            clustNumberBox.String = 'ID#';
        end
        
        % Activate action buttons
        activateActionButtons(true);
        
        myGlobalLoadInterrupt = 0; % Release lock
    end
    
    % function loadClustButtonCallback
    %   Additional level when requesting an electrode through a cluster ID
    function loadClustButtonCallback(source,callbackdata)
        c = str2num(clustNumberBox.String);
        if numel(c) ~= 1 || isnan(c) % Parsing
            statusBar.String = sprintf('Requested cluster %s not a number.',clustNumberBox.String);
            return;
        end
        % Get electrode number (and ID validity)
        el = backEndHandle.checkID(c);
        if el ~= -1 % (Invalid ID tag)
            elNumberBox.String = el-1;
            loadButtonCallback(source,callbackdata);
        end
    end
    
    % function refreshView()
    %   Refreshes all displays in the GUI
    function refreshView()
        makeColors(); % Generate colors for the cluster distribution
        clustMgmt.Data = getClusterData(); % Seek and update cluster data table
        refreshInfoRow(); % General information
        refreshGraphics(); % PC, rate, ACF plots
        refreshLowLeftPanels(); % EI and STA panels
        autoScalePCPlots(); % Autoscale of PC plots
    end
    
    % function makeColors
    %   Generates a set of colors for the loaded clusters
    %   Repartition is along an HSV H-circle
    %   With merged clusters closer in hue than others
    %
    %   Colors are stored in globalized clusterColors.RGB and .HSV
    function makeColors()
        % make HSV colors - depending on merge statuses (smaller steps)
        tmp = backEndHandle.statusRaw;
        idx = 1:backEndHandle.nClusters;
        tmp(tmp(:,1) ~= 2,2) = idx(tmp(:,1) ~= 2);
        tmp(tmp(:,1) == 2,2) = arrayfun(@(x) find(backEndHandle.displayIDs == x), tmp(tmp(:,1) == 2,2));
        tmp = tmp(:,2) - min(tmp(:,2)) + 1;
        [tmp,i] = sort(tmp); [~,j] = sort(i);
        tmp = tmp + 0.6 * (0:(backEndHandle.nClusters-1))';
        % Number    ----  is the merge color spacing coefficient. Increase or decrease
        % to have merge clusters be more or less close in color.
        tmp = tmp(j);
        clusterColors.HSV = [tmp ./ max(tmp),...
            ones(backEndHandle.nClusters,1)*[0.8 0.7]];
        if backEndHandle.nClusters > 0
            clusterColors.RGB = hsv2rgb(clusterColors.HSV); % convert to RGB
        else
            clusterColors.RGB = zeros(0,3);
        end
    end
    
    % function refreshGraphics
    %   Resamples and refreshes the contents of PCs, ACF and spike rate plots
    %   Is called from refreshGraphics and also the callback of the "reScatter" button
    function refreshGraphics(varargin)
        % varargin is (source,callbackdata) as used for button callback
        
        % Initialize cell arrays for plot handles
        PC123Plots = cell(backEndHandle.nClusters,1);
        PC45Plots = cell(backEndHandle.nClusters,1);
        ratePlots = cell(backEndHandle.nClusters,1);
        ACFPlots = cell(backEndHandle.nClusters,1);
        
        % Clear axes, put the holds to on for plot stacking
        cla(plot3D); hold(plot3D,'on');
        cla(PC45Plot); hold(PC45Plot,'on');
        cla(ratePlot); hold(ratePlot,'on');
        cla(ACFPlot); hold(ACFPlot,'on');
        
        % Display fraction - 1 if less than displayPoints points
        totSpikesForEl = sum(cellfun(@numel,backEndHandle.spikeTrains,'uni',true));
        displayFraction = min(1,displayPoints ./ totSpikesForEl);
            
        % Loop on clusters - plotting
        for c = 1:backEndHandle.nClusters
            % Compute the number of points for the cluster
            % Cluster relative fractions are preserved by the display sampling
            if displayFraction == 1 % use all points
                subsampleIndex = 1:size(backEndHandle.prjTrains{c},1);
            else % Randomly sample points
                subsampleIndex = randsample(size(backEndHandle.prjTrains{c},1),...
                    floor(displayFraction.*size(backEndHandle.prjTrains{c},1)));
            end
            
            % Generate plot colors point by point - Add a little shading along a given direction
            colorPlot = repmat(clusterColors.HSV(c,:),numel(subsampleIndex),1);
            colorPlot(:,2) = colorPlot(:,2) + 0.2*backEndHandle.prjTrains{c}(subsampleIndex,1)./...
                max(abs(backEndHandle.prjTrains{c}(subsampleIndex,1)));
            colorPlot = hsv2rgb(colorPlot);
            
            % Plot the 3D scatter
            PC123Plots{c} = scatter3(plot3D,...
                backEndHandle.prjTrains{c}(subsampleIndex,PC(1)),...
                backEndHandle.prjTrains{c}(subsampleIndex,PC(2)),...
                backEndHandle.prjTrains{c}(subsampleIndex,PC(3)),...
                9,...
                colorPlot,...
                'Visible',bool2onoff(clustMgmt.Data{c,3}));
            % Plot the 2D scatter
            PC45Plots{c} = scatter(PC45Plot,...
                backEndHandle.prjTrains{c}(subsampleIndex,PC(4)),...
                backEndHandle.prjTrains{c}(subsampleIndex,PC(5)),...
                9,...
                colorPlot,...
                'Visible',bool2onoff(clustMgmt.Data{c,3}));
            % Plot the spikeRate
            nBins = 100;
            [bars,centers] = hist(backEndHandle.spikeTrains{c} ./ 20000,nBins);
            ratePlots{c} = plot(ratePlot,...
                centers,bars * nBins * 20000/ backEndHandle.nSamples,...
                'color',clusterColors.RGB(c,:),...
                'Visible',bool2onoff(clustMgmt.Data{c,3}),...
                'LineWidth',1.5);
            % Plot the ACF - java computed
            correlator = edu.ucsc.neurobiology.vision.analysis.AutocorrelationCalculator.calculate(...
                backEndHandle.spikeTrains{c},100,0.5); % 200 msec extent @1msec resolution
            corrFun = correlator.toArray();
            corrFun = filter(0.125*ones(1,8),1,corrFun ./ sqrt(sum(corrFun.^2)));
            ACFPlots{c} = plot(ACFPlot,...
                correlator.getXValues,corrFun,...
                'color',clusterColors.RGB(c,:),...
                'Visible',bool2onoff(clustMgmt.Data{c,3}),...
                'LineWidth',1.5);
        end % for c = 1:backEndhandle.nClusters
    end % function refreshGraphics
    
    % function switchPC
    %   Changes the PC dimension displayed on the 3D and 2D PC plots
    %   callback to the switchPCOrderButton button
    function switchPC(~,~)
        params = inputdlg({'PCs for 3D plot ([X,Y,Z]):','PCs for 2D plot ([X,Y]):'},...
            'Input',[1 35],{prettyPrint(PC(1:3)),prettyPrint(PC(4:5))});
        try
            x = [str2num(params{1}),str2num(params{2})];
             catchvalidateattributes(x,{'numeric'},{'row','size',[1 5],'integer','>=',1,'<=',size(backEndHandle.prjTrains{1},2)},'','');
       
            statusBar.String = 'Invalid input - PC axes change aborted.';
            return;
        end
        PC = x;
        plot3D.Title.String = sprintf('Principal Components %u-%u-%u',PC(1),PC(2),PC(3));
        plot3D.XLabel.String = sprintf('PC %u',PC(1));
        plot3D.YLabel.String = sprintf('PC %u',PC(2));
        plot3D.ZLabel.String = sprintf('PC %u',PC(3));
        PC45Plot.Title.String = sprintf('Principal Components %u-%u',PC(4),PC(5));
        PC45Plot.XLabel.String = sprintf('PC %u',PC(4));
        PC45Plot.YLabel.String = sprintf('PC %u',PC(5));
        if backEndHandle.isDataReady;
            refreshView();
        end
        statusBar.String = 'PC Axes changed.';
    end
    
    % function refreshLowLeftPanels
    %   refreshes the lower-left quarter
    %   with EI, make-up EI, STA, EI distance matrix information
    %   Switches depending on single or multiple selection
    %   And depending if there exists EIs/STA for selection
    %
    % Single neuron selected:
    %   - neurons has EI:       - left shows EI
    %               else:       - left shows "No EI"
    %   - neuron has STA:       - right shows STA
    %               else:       - right shows "No STA"
    %
    % Multiple neurons selected
    %   - all selected have EI: - left shows make-up EI
    %                     else: - left shows "No EI"
    %   - right always shows EI distance matrix for the electrode
    function refreshLowLeftPanels()
        if numel(eiJPanel) > 0 % Clear previous ei java panel
            delete(eiJPanel); eiJPanel = [];
            delete(eiHPanel); eiHPanel = [];
        end
        if numel(staJPanel) > 0 % Clear previous sta java panel
            delete(staJPanel); staJPanel = [];
            delete(staHPanel); staHPanel = [];
        end
        % Get the cluster selection
        [k,c] = nClustSelected();
        
        % Single neuron selected - display EI and STA
        if k == 1
            leftColumns.Sizes(3) = -1; % Set EI/STA area size to auto
            bottomLeftAxes.Visible = 'off'; % Hide (unused) left axes
            bottomRightAxesBox.Visible = 'off'; % Hide the wrapper of EI dist matrix and its buttons
            
            % EI Panel side
            if numel(backEndHandle.eiFile) > 0 && numel(backEndHandle.eisLoaded{c}) > 0 % Single neuron has EI
                eiInfoString.String = 'Real neuron EI';
                [eiJPanel,eiHPanel] = ...
                    javacomponent(edu.ucsc.neurobiology.vision.neuronviewer.PhysiologicalImagePanel(...
                    backEndHandle.eisLoaded{c},...
                    [],2,...
                    backEndHandle.electrodeMap,...
                    backEndHandle.elLoaded - 1,...
                    '','',false,false,1.,true,false,true,40,true,3,0));
                eiHPanel.Parent = eiPanel;
                eiHPanel.Units = 'norm';
                eiHPanel.Position = [0 0 1 1];
                isEI = true;
            else % Single neuron doesn't have an EI
                isEI = false;
                eiInfoString.String = 'No EI for this neuron';
            end
            
            % STA Panel side
            if numel(backEndHandle.staFile) > 0 && numel(backEndHandle.stasLoaded{c}) > 0 % Single neuron has STA
                [staJPanel,staHPanel] = ...
                    javacomponent(edu.ucsc.neurobiology.vision.neuronviewer.STAPlotMaker.makeSTAPanel(...
                    backEndHandle.stasLoaded{c},...
                    false,1,1.0,false,true,false,...
                    backEndHandle.globalsFile,...
                    ''));
                staHPanel.Parent = staPanel;
                staHPanel.Units = 'norm';
                staHPanel.Position = [0 0 1 1];
                isSTA = true;
            else % Single neuron doesn't have an STA
                isSTA = false;
                bottomLeftAxes.Visible = 'off';
                noSTAString.Visible = 'on';
            end
            
            % Neither EI nor STA, fold the EI-STA panel to a single-liner
            if ~isEI && ~isSTA
                leftColumns.Sizes(3) = 26;
            end
            
        else % Multiple neurons selected - display EI distance and xCorr distance
            leftColumns.Sizes(3) = -1; % Expand EI-STA panel
            
            % EI Panel (Left)
            if k > 0 && numel(backEndHandle.eiFile) > 0 && all(cellfun(@(x) numel(x),backEndHandle.eisLoaded(c))) % All neurons have an EI
                eiInfoString.String = 'Make-up EI of selection';
                compoundEI = zeros(size(backEndHandle.eisLoaded{c(1)}));
                for neur = 1:numel(c)
                    compoundEI(1,:,:) = compoundEI(1,:,:) + backEndHandle.spikeCounts(c(neur)) .* backEndHandle.eisLoaded{c(neur)}(1,:,:);
                end
                compoundEI(1,:,:) = compoundEI(1,:,:) ./ sum(backEndHandle.spikeCounts(c));
                [eiJPanel,eiHPanel] = ...
                    javacomponent(edu.ucsc.neurobiology.vision.neuronviewer.PhysiologicalImagePanel(...
                    compoundEI,...
                    [],2,...
                    backEndHandle.electrodeMap,...
                    backEndHandle.elLoaded - 1,...
                    '','',false,false,1.,true,false,true,40,true,3,0));
                eiHPanel.Parent = eiPanel;
                eiHPanel.Units = 'norm';
                eiHPanel.Position = [0 0 1 1];
            else % Some selected neurons don't have an EI
                eiInfoString.String = 'No EI for this selection';
            end
            
            % Right side: display EI distance matrix
            noSTAString.Visible = 'off'; % hide the STA info string
            
            % Instantiate display in the axes and adjust all properties
            eiMatrix = imagesc(backEndHandle.EIdistMatrix(currentPermutation,currentPermutation),...
                'AlphaData',~isnan(backEndHandle.EIdistMatrix(currentPermutation,currentPermutation)),...
                'Parent',bottomRightAxes);
            bottomRightAxes.DataAspectRatio = [1 1 1];
            bottomRightAxes.XAxisLocation = 'top';
            bottomRightAxes.XTick = 1:backEndHandle.nClusters;
            bottomRightAxes.YTick = 1:backEndHandle.nClusters;
            bottomRightAxes.XTickLabel = num2cell(backEndHandle.displayIDs(currentPermutation));
            bottomRightAxes.YTickLabel = num2cell(backEndHandle.displayIDs(currentPermutation));
            bottomRightAxes.XTickLabelRotation = 35;
            bottomRightAxes.YTickLabelRotation = 35;
            bottomRightAxes.TickLength = [0,0];
            colorbar(bottomRightAxes,'eastoutside');
            colormap(bottomRightAxes,flipud(jet));
            bottomRightAxes.Position = [0.1    0    0.68    0.9];
            bottomRightAxes.CLim = [0,1];
            bottomRightAxesBox.Visible = 'on'; % Show the EI distance matrix
            
            % Plot a grid on the display
            nc = backEndHandle.nClusters;
            hold(bottomRightAxes,'on');
            plot(bottomRightAxes,repmat([0.5,nc + 0.5],nc-1,1)',((0.5+(1:(nc-1)))'*[1,1])','k-','linewidth',1);
            plot(bottomRightAxes,((0.5+(1:(nc-1)))'*[1,1])',repmat([0.5,nc + 0.5],nc-1,1)','k-','linewidth',1);
            hold(bottomRightAxes,'off');
        end
    end
    
    % function getClusterData
    % Seeks information from the backend
    % To fill the cluster data table
    function d = getClusterData()
        % HTML trick to show colored boxes in the table
        colorgen = @(color,text) sprintf('<html><table border=0 width=30 bgcolor=#%02X%02X%02X><TR><TD>%s</TD></TR> </table></html>',...
            ceil(255*color(1)),ceil(255*color(2)),ceil(255*color(3)),text);
        
        % Cluster IDs
        IDs = num2cell(backEndHandle.displayIDs);
        % Cluster colors
        colors = cellfun(@(x) colorgen(x,''),...
            mat2cell(clusterColors.RGB,ones(backEndHandle.nClusters,1)),...
            'uni',false);
        
        % Cluster statuses - build strings from the backend raw info
        status = cell(backEndHandle.nClusters,1);
        for c = 1:backEndHandle.nClusters
            switch backEndHandle.statusRaw(c,1)
                case -2
                    status{c} = 'Discard';
                case -1
                    status{c} = 'Unknown';
                case 0
                    status{c} = 'Keep';
                case 1
                    status{c} = 'Contam / Low count';
                case 2
                    status{c} = sprintf('Merge with %u',backEndHandle.statusRaw(c,2));
                case 3
                    [e,~] = backEndHandle.getElClust(backEndHandle.statusRaw(c,2));
                    status{c} = sprintf('Dup. of ID %i, El %u',backEndHandle.statusRaw(c,2),e-1);
            end
        end
        
        % Default display - show all except contaminated (1), low count (1) and discarded without
        % details (-2)
        display = num2cell(~or(backEndHandle.statusRaw(:,1) == 1,backEndHandle.statusRaw(:,1) == -2));
        
        contam = num2cell(backEndHandle.contaminationValues); % Contamination values
        spcount = num2cell(backEndHandle.spikeCounts); % Spike counts
        sprate = num2cell(backEndHandle.spikeCounts * 20000 / backEndHandle.nSamples); % Average spike rate
        comment = backEndHandle.comment; % Comment string (parsed from classification file)
        currentPermutation = 1:backEndHandle.nClusters; % reset EI dist matrix display permutation 
        
        % Define the data
        % ID - Color - Display - Status - Contam
        d = [IDs, colors, display, status, spcount, sprate, contam, comment];
    end
    
    % function selectorCallback
    %   Called when the show all/none buttons are clicked
    %   builds a callback data information for generalized handling in tableEditCallback
    function selectorCallback(source,callbackdata)
        if source == allButton % generate true NewData for "select all"
            fakeCallbackData.NewData = true;
        else % and false NewData for "select none"
            fakeCallbackData.NewData = false;
        end
        fakeCallbackData.Indices = [0,3]; % Indices(2) has to be 3 (selector column) for tableEditCallback
        
        % For each cluster except the last, unselect the checkbox, put the correct index, and call tableEditCallback
        % with FALSE opt argument (no display refresh)
        for c = 1:(backEndHandle.nClusters-1)
            clustMgmt.Data{c,3} = fakeCallbackData.NewData;
            fakeCallbackData.Indices(1) = c;
            tableEditCallback(source,fakeCallbackData,false);
        end
        
        % last cluster - call tableEditCallback WITH refresh option
        c = backEndHandle.nClusters;
        clustMgmt.Data{c,3} = fakeCallbackData.NewData;
        fakeCallbackData.Indices(1) = c;
        tableEditCallback(source,fakeCallbackData,true);
    end
    
    % function tableEditCallback
    %   callback of cluster data table
    %   called when a cluster display checkbox is switched
    %   shows/hides corresponding plots
    %   If we go from/to multiple selection to single selection, also
    %   refreshes lower left area accordingly
    %
    % Takes-in optional third boolean argument - to inhibit (false) or force (true)
    % display refreshing, to prevent useless operations when all clusters are
    % selected or unselected at once by the corresponding buttons
    function tableEditCallback(source,callbackdata,varargin)
        if ~(callbackdata.Indices(2) == 3) % Only third column can be a valid call
            throw(MException('','ClusterEditGUI:tableEditCallback - Call from not a tickbox'));
        end
        % Show-hide all 4 by-cluster plots, and put activated cluster on top of stack
        PC123Plots{callbackdata.Indices(1)}.Visible = bool2onoff(callbackdata.NewData);
        uistack(PC123Plots{callbackdata.Indices(1)},'top');
        PC45Plots{callbackdata.Indices(1)}.Visible = bool2onoff(callbackdata.NewData);
        uistack(PC45Plots{callbackdata.Indices(1)},'top');
        ratePlots{callbackdata.Indices(1)}.Visible = bool2onoff(callbackdata.NewData);
        uistack(ratePlots{callbackdata.Indices(1)},'top');
        ACFPlots{callbackdata.Indices(1)}.Visible = bool2onoff(callbackdata.NewData);
        uistack(ACFPlots{callbackdata.Indices(1)},'top');
        
        % Filter if caller requests manual refresh or not - Allows to skip refresh in case
        % of successive calls from "(Un)Select All" buttons.
        if nargin == 2 || (nargin == 3 && varargin{1} == true)
            refreshLowLeftPanels();
        end
    end
    
    % function autoScalePCPlots
    %   auto scales the axes of the PC plots to their currently visible contents
    %   callback of the "rescale axes" button
    %   Otherwise, these plots hold their current axes
    function autoScalePCPlots(varargin)
        listObj = [plot3D;PC45Plot];
        % The drawnow MUST BE UNINTERRUPTED
        % Otherwise (or if there is no drawnow)
        % This very callback interrupts the drawnow and sets back lim modes to
        % manual. When the drawnow resumes, LimModes are manual so there is
        % nothing to be refreshed...
        set(listObj,'Interruptible','off');
        set(listObj,'XLimMode','auto');
        set(listObj,'YLimMode','auto');
        set(listObj,'ZLimMode','auto');
        drawnow;
        set(listObj,'XLimMode','manual');
        set(listObj,'YLimMode','manual');
        set(listObj,'ZLimMode','manual');
        set(listObj,'Interruptible','on');
    end
    
    % function nClustSelected()
    %   returns the number of clusters currently selected in the cluster data table
    %   and the corresponding selection
    function [k,c] = nClustSelected()
        k = sum(cell2mat(clustMgmt.Data(:,3)));
        c = find(cell2mat(clustMgmt.Data(:,3)));
    end
    
    % function refreshInfoRow()
    %   refreshes the information of the global information row
    %   with what is currently loaded in the backend 
    function refreshInfoRow()
        infoRow.String = sprintf('El %u:  %u Total Clusters  |  %u Contam/Low Count  |  %u Locally Merged  |  %u Global Duplicates  |  %u Final',...
            backEndHandle.elLoaded-1,...
            backEndHandle.nClusters,...
            nnz(backEndHandle.statusRaw(:,1) == 1),...
            nnz(backEndHandle.statusRaw(:,1) == 2),...
            nnz(backEndHandle.statusRaw(:,1) == 3),...
            nnz(backEndHandle.statusRaw(:,1) == 0));
    end
    
    % function defaultPermutationCallback(...)
    %   callback of the "Default" button of the EI distance matrix display
    %   sets the display permutation ot natural ordering 
    function defaultPermutationCallback(varargin)
        currentPermutation = 1:backEndHandle.nClusters;
        refreshLowLeftPanels();
    end
    
    % function optimalPermutationCallback(...)
    %   callback of the "Optimal" button of the EI distance matrix display
    %   sets the display permutation to the optimal one with the metric
    %   defined in optimalBlockDiagPerm()
    %   Discard the clusters with missing EIs from the display before performing this.
    function optimalPermutationCallback(varargin)
        validCol = find(cellfun(@(x) numel(x) > 0, backEndHandle.eisLoaded));
        currentPermutation = validCol(optimalBlockDiagPerm(0.5 * (2-backEndHandle.EIdistMatrix(validCol,validCol))));
        refreshLowLeftPanels();
    end
    
    % function customPermutationCallback(...)
    %   callback of the "Custom" button of the EI distance matrix display
    %   Opens a prompt window for the user to provide a custom list of IDs to display in that order
    function customPermutationCallback(varargin)
        s = prettyPrint(backEndHandle.displayIDs(currentPermutation));
        request = inputdlg('Enter cluster ordering:','Input',[1 100],{s});
        if numel(request) == 0
            return % User clicked cancel
        end
        try % try to evaluate user input as valid matlab content, with valid ID number for current electrode
            eval(['currentPermutation = ',request{1},';']);
            currentPermutation = arrayfun(@(x) find(backEndHandle.displayIDs == x), currentPermutation);
            statusBar.String = 'EI display permutation changed';
            refreshLowLeftPanels();
        catch
            statusBar.String = 'Invalid permutation entered';
        end
    end
    
    % function clearEditsCallback
    % Removes all edits being applied/previewed on the currently displayed electrode
    % and reloads electrode data from scratch.
    % Callback of the clear edits button
    function clearEditsCallback(source,~)
        if source == clearEditActions
            editHandler.clearSession();
        elseif source == deleteEditActions
            editHandler.clearLoad();
            editHandler.requestClearHard();
        end
        backEndHandle.reloadSameData();
        for i = 1:size(editHandler.editList,1)
            backEndHandle.localApplyAction(editHandler.editList{i,1},editHandler.editList{i,2},editHandler.editList{i,3});
        end
        if numel(editHandler.editList) > 0
            editHandler.openWindow();
        end
        refreshView();
    end
    
    % function saveEditsCallback
    % prompts user for validation
    % save edits to .edit.mat file
    function saveEditsCallback(varargin)
        editHandler.saveEdits();
    end
    
    % function closeWinCallback
    %   is the DeleteFcn of the principal window
    %   will prompt for saving pending edits then close the editHandler
    %   window in case it's open
    function closeWinCallback(varargin)
        saveEditsCallback();
        editHandler.closeWindow();
    end
    
    % function openEditHandlerCallback
    %   opens the editHandler window
    %   callback to corresponding button
    function openEditHandlerCallback(~,~)
        editHandler.openWindow();
    end
    
    % function editCallback
    %   callback of all edition buttons
    %   Checks the action nature and applies appropriate behavior
    %   updates the editHandles and calls the backend for applying changes
    function editCallback(~,~, action)
        [k,c] = nClustSelected();
        c = backEndHandle.displayIDs(c);
        s = prettyPrint(c);
        try
            exc = [];
            switch action
                case EditAction.KEEP
                    params = inputdlg({'List of IDs to elevate:'},'Input',[1 70],{s});
                    params{1} = str2num(params{1});
                case EditAction.DELETE
                    params = inputdlg({'List of IDs to delete:'},'Input',[1 70],{s});
                    params{1} = str2num(params{1});
                case EditAction.MERGE
                    params = inputdlg({'List of IDs to merge together:'},'Input',[1 70],{s});
                    params{1} = str2num(params{1});
                case EditAction.RECLUSTER
                    params = inputdlg({'List of IDs to recluster:','Number of clusters (0 for automatic):','Configuration tag (empty for default):'},...
                        'Input',[1 70],{s,num2str(k),'LIGHT_CLUST'});
                    params{1} = str2num(params{1});
                    params{2} = str2double(params{2});
                case EditAction.SHRINK
                    params = inputdlg({'List of IDs to shrink:','Target fraction of remaining spikes (1 to skip):',...
                        'Target contamination (inf to skip):'},'Input',[1 70],{s,'1','inf'});
                    params{1} = str2num(params{1});
                    params{2} = str2double(params{2});
                    params{3} = str2double(params{3});
                otherwise
                    % report unhandled case out of try/catch block
                    exc = MException('','ClusterEditGUI:editCallback - Unhandled EditAction in switch statement.');
            end
        catch
            params = {[]};
        end
        if numel(exc) > 0
            throw(exc)
        end
        if numel(params) == 0
            return % User clicked cancel
        end
        [v,m] = action.checkParameters(params);
        if v == 0
            statusBar.String = 'Invalid action parameters. See console for details';
            fprintf([m,'\n']);
            return;
        end
        
        statusBar.String = [action.char,' request acknowledged, computing...']; drawnow;
        % First check backend for validity of execution
        % Then add to edithandler and open it
        [s,data] = backEndHandle.localApplyAction(action,params);
        if s ~= 0 % Invalid parameters regarding the backend data, aborting.
            statusBar.String = data;
            return;
        end
        
        editHandler.addAction(action,params,data,true);
        editHandler.openWindow();
        
        statusBar.String = 'Computation done. Refreshing view...'; drawnow;
        refreshView();
        statusBar.String = [action.char,' edition done.']; drawnow;
        
    end
    
    % function activateActionButtons
    %   activate or deactivate the edit action buttons
    %   Input:
    %       onOrOff: boolean, true to activate, false to deactivate
    function activateActionButtons(onOrOff)
        if backEndHandle.typeIsSpectra % Edition not allowed in vision mode.
            set(openEditHandler,'Enable',bool2onoff(onOrOff));
            set(clearEditActions,'Enable',bool2onoff(onOrOff));
            set(deleteEditActions,'Enable',bool2onoff(onOrOff));
            set(saveEditActions,'Enable',bool2onoff(onOrOff));
            set(undoLastEditButton,'Enable',bool2onoff(onOrOff));
            cellfun(@(x) set(x,'Enable',bool2onoff(onOrOff)), actionButtons);
        end
    end
    
    % function undoLastEdit
    %   undo the last edit action
    %   callback of the undoLastEditButton pushbutton
    function undoLastEdit(~,~)
        if size(editHandler.editList,1) == 0
            statusBar.String = 'No action to remove';
            return;
        end
        if ~editHandler.unsavedActions
            statusBar.String = 'Cannot Remove saved actions this way. Please use -DEL ALL- then -SAVE-.';
            return;
        end
        editHandler.removeLastAction();
        statusBar.String = 'Canceled latest action.';
        backEndHandle.reloadSameData();
        for i = 1:size(editHandler.editList,1)
            backEndHandle.localApplyAction(editHandler.editList{i,1}, editHandler.editList{i,2}, editHandler.editList{i,3});
        end
        refreshView();
    end
    
    %% FINISH GENERATING THE GUI %%
    if backEndHandle.isDataReady;
        refreshView();
        elNumberBox.String = num2str(backEndHandle.elLoaded);
    else
        activateActionButtons(false);
    end
    frontEndHandle.Visible = 'on';
    varargout{1} = frontEndHandle;
    varargout{2} = backEndHandle;
end

%% ---------------------------------------------------------- %%
% MISC UTILITY FUNCTIONS
%---------------------------------------------------------------

% function bool2onoff()
% Takes a boolean and returns appropriate 'on' on 'off' string
function s = bool2onoff(b)
    if b == 0
        s = 'off';
    elseif b == 1
        s = 'on';
    else
        throw(MException('','bool2onoff: not 0 or 1'));
    end
end

% function bool2tf()
% Takes a boolean and returns appropriate 'true' or 'false' string
function s = bool2tf(b)
    if b == 0
        s = 'false';
    elseif b == 1
        s = 'true';
    else
        throw(MException('','bool2tf: not 0 or 1'));
    end
end

% function customizeTools()
%   removes all the icons we don't want in the default figure toolbar
%   We keep drag, zoom in/out, data points and rotate.
function [remainingTools,toolbarHandle] = customizeTools(figureHandle)
        allTools = findall(figureHandle);
        deleteTags = {'Plottools*',...
            'Annotation*',...
            'DataManager*','Standard*'};
        for i = 1:numel(deleteTags);
            t = findobj(allTools,'-regexp','Tag',deleteTags{i});
            set(t,'Separator','off');
            set(t,'Visible','off');
        end
        t = findobj(allTools,'Tag','Exploration.ZoomIn');
        t.Separator = 'off';
        t = findobj(allTools,'TooltipString','Brush/Select Data');
        set(t,'Visible','off');
        
        remainingTools = findobj(allTools,'-regexp','Tag','Exploration*','-and','Visible','on');
        toolbarHandle = findobj(allTools,'Tag','FigureToolBar');
end