function [backEndHandle,frontEndHandle] = ClusterEditGUI(datasetFolder,varargin)
    % function ClusterEditGUI
    % Some comments here
    
    %% Instantiation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    narginchk(1,2); % Existing back-end may be passed in varargin{1}
    if nargin == 1
        backEndHandle = ClusterEditBackend(datasetFolder);
    else
        backEndHandle = varargin{1};
    end
    % Front-end main figure
    frontEndHandle = figure( 'Name', 'Cluster Editor 0.0', ...
        'MenuBar', 'none', ...
        'Toolbar', 'figure', ...
        'NumberTitle', 'off',...
        'Visible', 'off',...
        'OuterPosition',get(groot,'Screensize') + [60 60 -120 -90]);
    customizeTools(frontEndHandle); % Remove useless buttons
    
    %% Some configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    spacerWidth = 10;
    displayPoints = 8000;
    clusterColors = zeros(0,3);
    
    %% GUI layout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main box and children
    mainLayout = uiextras.HBoxFlex(...
        'Parent',frontEndHandle,...
        'Position',[0 0 1 1],...
        'Spacing',2);
    % Left VBox - Menu, title, table, EI/STA displays
    leftColumns = uiextras.VBox(...
        'Parent',mainLayout,...
        'Spacing',spacerWidth);
    % Right VBox - PC, ACF and spRate
    graphLayout = uiextras.VBoxFlex(...
        'Parent',mainLayout,...
        'Spacing',spacerWidth);
    mainLayout.Sizes = [-1 -1];
    % Graph 3D + view buttons strip
    graph3DBox = uiextras.HBox(...
        'Parent',graphLayout,...
        'Spacing',spacerWidth);
    % View button column
    %{
    % Removing button column - available w/ right click while rotate tool
    buttonColumn = uiextras.VBox(...
        'Parent',graph3DBox,...
        'Spacing',spacerWidth);
    buttonDefaultView = uicontrol(...
        'Parent',buttonColumn,...
        'Style', 'pushbutton',...
        'String','Default',...
        'fontsize',11,...
        'callback',@view3DCallback);
    buttonXY = uicontrol(...
        'Parent',buttonColumn,...
        'Style', 'pushbutton',...
        'String','X-Y',...
        'fontsize',11,...
        'callback',@view3DCallback);
    buttonXZ = uicontrol(...
        'Parent',buttonColumn,...
        'Style', 'pushbutton',...
        'String','X-Z',...
        'fontsize',11,...
        'callback',@view3DCallback);
    buttonYZ = uicontrol(...
        'Parent',buttonColumn,...
        'Style', 'pushbutton',...
        'String','Y-Z',...
        'fontsize',11,...
        'callback',@view3DCallback);
    buttonColumn.Sizes = [24 24 24 24];
    %}
    % Panel for 3D box axes
    supportPanel3D = uipanel('Parent',graph3DBox);
    PC123Plots = {}; % handle to subplots handle defined in load callback
    plot3D = axes(...
        'Parent',supportPanel3D,...
        'ClippingStyle','rectangle',...
        'View',[45,15],...
        'XGrid','on','YGrid','on','ZGrid','on',...
        'DataAspectRatio',[1 1 1],...
        'XLimMode','manual','YLimMode','manual','ZLimMode','manual');
    plot3D.Title.String = 'Principal Components 1-2-3';
    plot3D.XLabel.String = 'PC 1';
    plot3D.YLabel.String = 'PC 2';
    plot3D.ZLabel.String = 'PC 3';
    axis(plot3D,'tight');
    graph3DBox.Sizes = -1; %[50 -1]; % Sizes with button column
    bottomGraphs = uiextras.HBoxFlex(...
        'Parent',graphLayout,...
        'Spacing',spacerWidth);
    graphLayout.Sizes = [-3 -2];
    bottomGraphsLeftBox = uiextras.VBox(...
        'Parent',bottomGraphs,...
        'Spacing',0);
    bottomGraphsRightPanel = uipanel(...
        'Parent',bottomGraphs);
    bottomGraphs.Sizes = [-1 -1];
    
    % Three 2D plots at bottom right
    % Panels
    ratePanel = uipanel(...
        'Parent',bottomGraphsLeftBox);
    ratePlots = {}; % global handle to subplots handle defined in load callback
    ratePlot = axes(...
        'Parent',ratePanel,...
        'XGrid','on','YGrid','on',...
        'OuterPosition',[0 0 1 1]);
    ratePlot.Title.String = 'Spike rates'; ratePlot.Title.FontSize = 9;
    ratePlot.XLabel.String = 'Time (sec)'; ratePlot.XLabel.FontSize = 9;
    ratePlot.YLabel.String = 'Spike rate (Hz)'; ratePlot.YLabel.FontSize = 9;
    ratePlot.XLim = [0, backEndHandle.nSamples / 20000 + 1];
    
    ACFPanel = uipanel(...
        'Parent',bottomGraphsLeftBox);
    ACFPlots = {}; % global handle to subplots handle defined in load callback
    ACFPlot = axes(...
        'Parent',ACFPanel,...
        'XGrid','on','YGrid','on',...
        'OuterPosition',[0 0 1 1]);
    ACFPlot.Title.String = 'ACF'; ACFPlot.Title.FontSize = 9;
    ACFPlot.XLabel.String = 'Time \Delta (msec)'; ACFPlot.XLabel.FontSize = 9;
    ACFPlot.YLabel.String = 'Autocorr. (pair fraction/msec \Delta)'; ACFPlot.YLabel.FontSize = 9;
    ACFPlot.XLim = [0, 101];
    bottomGraphsLeftBox.Sizes = [-1 -1];
    
    % 4th col bottom
    PC45Plots = {}; % global handle to subplots handle defined in load callback
    PC45Plot = axes(...
        'Parent',bottomGraphsRightPanel,...
        'XGrid','on','YGrid','on',...
        'OuterPosition',[0 0 1 1],...
        'XLimMode','manual','YLimMode','manual','ZLimMode','manual');
    PC45Plot.Title.String = 'Principal Components 4-5';
    PC45Plot.XLabel.String = 'PC 4';
    PC45Plot.YLabel.String = 'PC 5';
    axis(PC45Plot,'tight');
    
    datasetName = uicontrol(...
        'Parent',leftColumns,...
        'Style','text',...
        'fontsize',12,...
        'String',datasetFolder);
    menu = uiextras.VBox(...
        'Parent',leftColumns,...
        'Spacing',0);
    eistaBox = uiextras.HBox('Parent',leftColumns,...
        'Spacing',0);
    statusBar = uicontrol(...
        'Parent',leftColumns,...
        'Style','text',...
        'fontsize',10,...
        'String','Welcome to Cluster Editor 0.0 !',...
        'HorizontalAlignment','left');
    
    % Pass status bar handle to backend
    backEndHandle.statusBarHandle = statusBar;
    
    leftColumns.Sizes = [24 -1 0 24];
    
    
    % Left column menu - described by rows.
    loadRow = uiextras.HBox(...
        'Parent',menu,...
        'Spacing',6,...
        'Padding',6);
    
    elNumberBox = uicontrol(...
        'Parent',loadRow,...
        'Style', 'edit',...
        'Max',1,'Min',1,...
        'String','EL#',...
        'fontsize',11,...
        'callback',@loadButtonCallback);
    
    ppmmButtonBox = uiextras.VBox(...
        'Parent',loadRow,...
        'Spacing',0,...
        'Padding',0);
    
    ppButton = uicontrol(...
        'Parent',ppmmButtonBox,...
        'Style', 'pushbutton',...
        'String',char(hex2dec('25B2')),...
        'fontsize',8,...
        'callback',@loadButtonCallback);
    
    mmButton = uicontrol(...
        'Parent',ppmmButtonBox,...
        'Style', 'pushbutton',...
        'String',char(hex2dec('25BC')),...
        'fontsize',8,...
        'callback',@loadButtonCallback);
    ppmmButtonBox.Sizes = [-1 -1];
    
    loadButton = uicontrol(...
        'Parent',loadRow,...
        'Style', 'pushbutton',...
        'String','Load',...
        'fontsize',11,...
        'callback',@loadButtonCallback);
    
    clustNumberBox = uicontrol(...
        'Parent',loadRow,...
        'Style', 'edit',...
        'Max',1,'Min',1,...
        'String','ID#',...
        'fontsize',11,...
        'callback',@loadClustButtonCallback);
    
    loadClusterButton = uicontrol(...
        'Parent',loadRow,...
        'Style', 'pushbutton',...
        'String','Load',...
        'fontsize',11,...
        'callback',@loadClustButtonCallback);
    
    fillerLoadRow = uicontrol(...
        'Parent',loadRow,...
        'Style','text',...
        'String','');
    
    % finish loadRow
    loadRow.Sizes = [34, 20, 50, 45, 60, -1];
    
    % selector row w/ "(un)select all"
    selectorRow = uiextras.HBox(...
        'Parent',menu,...
        'Spacing',6,...
        'Padding',6);
    
    allButton = uicontrol(...
        'Parent',selectorRow,...
        'Style', 'pushbutton',...
        'String','Select All',...
        'fontsize',11,...
        'callback',@selectorCallback);
    noneButton = uicontrol(...
        'Parent',selectorRow,...
        'Style', 'pushbutton',...
        'String','Unselect All',...
        'fontsize',11,...
        'callback',@selectorCallback);
    reScatterButton = uicontrol(...
        'Parent',selectorRow,...
        'Style', 'pushbutton',...
        'String','Re-scatter',...
        'fontsize',11,...
        'callback',@refreshGraphics);
    autoscaleButton = uicontrol(...
        'Parent',selectorRow,...
        'Style', 'pushbutton',...
        'String','Rescale Axes',...
        'fontsize',11,...
        'callback',@autoScalePCPlots);
    fillerSelectorRow = uicontrol(...
        'Parent',selectorRow,...
        'Style','text',...
        'String','');
    
    selectorRow.Sizes = [75 90 75 95 -1];
    % finish selectorRow
    
    % Column names and column format
    columnName = {'ID','','Disp','Status','#Spikes','Rate (Hz)','Contam','Classification'};
    columnFormat = {'numeric','char','logical','char','numeric','numeric','numeric','char'};
    columnEdit = [false, false, true, false, false, false, false, false];
    px = 10;colWidth =  {px*4, px*4, px*3, px*12, px*5, px*6, px*6, px*10};
    d = cell(0,8);
    clustMgmt = uitable(...
        'Parent',menu,...
        'ColumnName', columnName,...
        'ColumnFormat', columnFormat,...
        'ColumnEditable', columnEdit,...
        'RowName',[],...
        'CellEditCallback',@tableEditCallback);
    jscrollpane = findjobj(clustMgmt);
    jtable = jscrollpane.getViewport.getView;
    
    % Now turn the JIDE sorting on
    jtable.setSortable(true);
    jtable.setAutoResort(true);
    jtable.setColumnAutoResizable(true);
    jtable.setMultiColumnSortable(true);
    jtable.setPreserveSelectionsAfterSorting(true);
    clustMgmt.ColumnWidth = colWidth;
    clustMgmt.Data = d;
    % Finish left menu
    menu.Sizes = [34 34 -1];
    
    
    eiBox = uiextras.VBox(...
        'Parent',eistaBox,...
        'Padding',0,'Spacing',0);
    noEIString = uicontrol(...
        'Parent',eiBox,...
        'Style','text',...
        'fontsize',12,...
        'String','EI Display',...
        'Visible','on');
    noEIString.Units = 'Norm';
    noEIString.Position = [0 0 1 1];
    eiPanel = uipanel('Parent',eiBox,...
        'BorderType','none');
    eiJPanel = []; eiHPanel = []; % Globalize the java panels
    bottomLeftAxes = axes('Parent',eiPanel,...
        'Visible','off');
    
    eiBox.Sizes = [26 -1];
    
    staPanel = uipanel('Parent',eistaBox,...
        'BorderType','none');
    staJPanel = []; staHPanel = []; % Globalize the java panels
    
    bottomRightAxes = axes('Parent',staPanel,...
        'Visible','off');
    noSTAString = uicontrol(...
        'Parent',staPanel,...
        'Style','text',...
        'fontsize',12,...
        'String','No STA for this neuron',...
        'Visible','off');
    noSTAString.Units = 'Norm';
    noSTAString.Position = [0 0 1 1];
    
    eistaBox.Sizes = [-1 -1];
    
    %---------------------------------------------------------------
    % CALLBACKS AND SUBFUNCTIONS
    %---------------------------------------------------------------
    myGlobalLoadInterrupt = 0;
    function loadButtonCallback(source,callbackdata)
        myGlobalLoadInterrupt = 1;
        
        e = str2num(elNumberBox.String) + 1; % String is in JAVA numbering, e is in matlab numbering.
        if numel(e) ~= 1 || isnan(e)
            statusBar.String = sprintf('Requested electrode %s not a number.',elNumberBox.String);
            return;
        end
        if source == ppButton
            elNumberBox.String = num2str(e);
            e = e+1;
        end
        if source == mmButton
            elNumberBox.String = num2str(e-2);
            e = e-1;
        end
        statusBar.String = sprintf('Loading electrode %s...',elNumberBox.String);
        % Grab the interruption queue if there is one.
        % If there is none, disable interruption while uploading backend
        pause(.5); % Give us time to get interrupted (avoid redudancy queuing in case of frantic clicking on ++/-- buttons)
        % If interrupted, do not resume
        if ~myGlobalLoadInterrupt
            return;
        end
        source.Interruptible = 'off';
        
        status = backEndHandle.loadEl(e);
        if status ~= 0
            return;
        end
        
        % Allow for interruption again - before refreshing front end
        % And grab interruption queue
        source.Interruptible = 'on';
        drawnow;
        if ~myGlobalLoadInterrupt
            return;
        end
        % If interrupted, do not resume
        
        statusBar.String = 'Refreshing View...'; drawnow;
        refreshView();
        statusBar.String = sprintf('Displaying electrode %u.',e-1);
        if source == loadButton || source == elNumberBox || source == ppButton || source == mmButton
            clustNumberBox.String = 'ID#';
        end
        myGlobalLoadInterrupt = 0;
    end
    
    function loadClustButtonCallback(source,callbackdata)
        statusBar.String = sprintf('Loading neuron ID %s...',clustNumberBox.String);
        c = str2num(clustNumberBox.String);
        if numel(c) ~= 1 || isnan(c)
            statusBar.String = sprintf('Requested cluster %s not a number.',clustNumberBox.String);
            return;
        end
        el = backEndHandle.checkID(c);
        if el ~= -1
            elNumberBox.String = el-1;
            loadButtonCallback(source,callbackdata);
        end
    end
    
    function refreshView()
        makeColors();
        clustMgmt.Data = getClusterData();
        refreshGraphics();
        refreshLowLeftPanels();
        autoScalePCPlots();
    end
    
    function makeColors()
        % make HSV colors
        %clusterColors.HSV = [randperm(backEndHandle.nClusters)' / backEndHandle.nClusters,...
        %    ones(backEndHandle.nClusters,1)*[0.8 0.7]];
        tmp = backEndHandle.statusRaw;
        tmp(tmp(:,1) ~= 2,2) = backEndHandle.displayIDs(tmp(:,1) ~= 2);
        tmp = tmp(:,2) - min(tmp(:,2)) + 1;
        [tmp,i] = sort(tmp); [~,j] = sort(i);
        tmp = tmp + 0.75 * (0:(backEndHandle.nClusters-1))';
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
    
    function refreshGraphics(varargin)
        % varargin to handle (source,callbackdata) as used for button callback
        % Go for plots
        PC123Plots = cell(backEndHandle.nClusters,1);
        PC45Plots = cell(backEndHandle.nClusters,1);
        ratePlots = cell(backEndHandle.nClusters,1);
        ACFPlots = cell(backEndHandle.nClusters,1);
        
        % Put the holds to off for the first draw, then put them back to on at end of loop
        cla(plot3D); hold(plot3D,'on');
        cla(PC45Plot); hold(PC45Plot,'on');
        cla(ratePlot); hold(ratePlot,'on');
        cla(ACFPlot); hold(ACFPlot,'on');
        for c = 1:backEndHandle.nClusters
            % prjTrain is subsampled to the number of viewpoints
            % set in the header of the back end class.
            % Cluster relative fractions are conserved.
            displayFraction = min(1,displayPoints ./ size(backEndHandle.prjLoaded,1));
            if displayFraction == 1
                subsampleIndex = 1:size(backEndHandle.prjTrains{c},1);
            else
                subsampleIndex = randsample(size(backEndHandle.prjTrains{c},1),...
                    floor(displayFraction.*size(backEndHandle.prjTrains{c},1)));
            end
            
            % TODO rework the HSV shading
            colorPlot = repmat(clusterColors.HSV(c,:),numel(subsampleIndex),1);
            colorPlot(:,2) = colorPlot(:,2) + 0.2*backEndHandle.prjTrains{c}(subsampleIndex,1)./...
                max(abs(backEndHandle.prjTrains{c}(subsampleIndex,1)));
            colorPlot = hsv2rgb(colorPlot);
            
            PC123Plots{c} = scatter3(plot3D,...
                backEndHandle.prjTrains{c}(subsampleIndex,1),...
                backEndHandle.prjTrains{c}(subsampleIndex,2),...
                backEndHandle.prjTrains{c}(subsampleIndex,3),...
                9,...
                colorPlot,...
                'Visible',bool2onoff(clustMgmt.Data{c,3}));
            PC45Plots{c} = scatter(PC45Plot,...
                backEndHandle.prjTrains{c}(subsampleIndex,4),...
                backEndHandle.prjTrains{c}(subsampleIndex,5),...
                9,...
                colorPlot,...
                'Visible',bool2onoff(clustMgmt.Data{c,3}));
            nBins = 100;
            [bars,centers] = hist(backEndHandle.spikeTrains{c} ./ 20000,nBins);
            ratePlots{c} = plot(ratePlot,...
                centers,bars * nBins * 20000/ backEndHandle.nSamples,...
                'color',clusterColors.RGB(c,:),...
                'Visible',bool2onoff(clustMgmt.Data{c,3}),...
                'LineWidth',1.5);
            correlator = edu.ucsc.neurobiology.vision.analysis.AutocorrelationCalculator.calculate(...
                backEndHandle.spikeTrains{c},100,0.5); % 200 msec extent @1msec resolution
            corrFun = correlator.toArray();
            corrFun = filter(0.125*ones(1,8),1,corrFun ./ sqrt(sum(corrFun.^2)));
            ACFPlots{c} = plot(ACFPlot,...
                correlator.getXValues,corrFun,...
                'color',clusterColors.RGB(c,:),...
                'Visible',bool2onoff(clustMgmt.Data{c,3}),...
                'LineWidth',1.5);
            
        end
    end
    
    % Refresh EI and STA display
    % Show if 1 neuron selected
    % Hide otherwise
    function refreshLowLeftPanels()
        
        if numel(eiJPanel) > 0 % Clear previous ei java panel
            delete(eiJPanel); eiJPanel = [];
            delete(eiHPanel); eiHPanel = [];
        end
        if numel(staJPanel) > 0 % Clear previous ei java panel
            delete(staJPanel); staJPanel = [];
            delete(staHPanel); staHPanel = [];
        end
        
        [k,c] = nClustSelected();
        if k == 1 % Single neuron selected - display EI and STA
            leftColumns.Sizes(3) = -1;
            bottomLeftAxes.Visible = 'off';
            bottomRightAxes.Visible = 'off';
            % EI Panel
            if numel(backEndHandle.eiFile) > 0 && numel(backEndHandle.eisLoaded{c}) > 0
                noEIString.String = 'Real neuron EI';
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
            else
                isEI = false;
                noEIString.String = 'No EI for this neuron';
            end
            
            % STA Panel
            if numel(backEndHandle.staFile) > 0 && numel(backEndHandle.stasLoaded{c}) > 0
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
            else
                isSTA = false;
                bottomLeftAxes.Visible = 'off';
                noSTAString.Visible = 'on';
            end
            
            if ~isEI && ~isSTA
                leftColumns.Sizes(3) = 26;
            end
            
        else % Multiple neurons selected - display EI distance and xCorr distance
            leftColumns.Sizes(3) = -1;
            % EI Panel (Left)
            
            if k > 0 && numel(backEndHandle.eiFile) > 0 && all(cellfun(@(x) numel(x),backEndHandle.eisLoaded(c)))
                noEIString.String = 'Make-up EI of selection';
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
            else
                noEIString.String = 'No EI for this selection';
            end
            
            % EI Dist Panel - Right
            noSTAString.Visible = 'off';
            
            eiMatrix = imagesc(backEndHandle.EIdistMatrix,...
                'AlphaData',~isnan(backEndHandle.EIdistMatrix),...
                'Parent',bottomRightAxes);
            title(bottomRightAxes,'Pairwise EI distance');
            bottomRightAxes.DataAspectRatio = [1 1 1];
            bottomRightAxes.XAxisLocation = 'top';
            bottomRightAxes.XTick = 1:backEndHandle.nClusters;
            bottomRightAxes.YTick = 1:backEndHandle.nClusters;
            bottomRightAxes.XTickLabel = num2cell(backEndHandle.displayIDs);
            bottomRightAxes.YTickLabel = num2cell(backEndHandle.displayIDs);
            bottomRightAxes.XTickLabelRotation = 35;
            bottomRightAxes.YTickLabelRotation = 35;
            bottomRightAxes.TickLength = [0,0];
            colorbar(bottomRightAxes,'eastoutside');
            bottomRightAxes.Position = [0.1    0    0.68    0.9];
            bottomRightAxes.CLim = [0,1];
            bottomRightAxes.Visible = 'on';
            
            nc = backEndHandle.nClusters;
            hold(bottomRightAxes,'on');
            plot(bottomRightAxes,repmat([0.5,nc + 0.5],nc-1,1)',((0.5+(1:(nc-1)))'*[1,1])','k-','linewidth',1);
            plot(bottomRightAxes,((0.5+(1:(nc-1)))'*[1,1])',repmat([0.5,nc + 0.5],nc-1,1)','k-','linewidth',1);
            hold(bottomRightAxes,'off');
        end
    end
    
    function d = getClusterData()
        % HTML trick
        colorgen = @(color,text) sprintf('<html><table border=0 width=30 bgcolor=#%02X%02X%02X><TR><TD>%s</TD></TR> </table></html>',...
            ceil(255*color(1)),ceil(255*color(2)),ceil(255*color(3)),text);
        
        %         % Column names and column format
        %         columnname = {'ID','','Disp','Status','Contam'};
        %         columnformat = {'char','char','logical','char','shortg'};
        IDs = num2cell(backEndHandle.displayIDs);
        colors = cellfun(@(x) colorgen(x,''),...
            mat2cell(clusterColors.RGB,ones(backEndHandle.nClusters,1)),...
            'uni',false);
        status = cell(backEndHandle.nClusters,1);
        % Parse neuron statuses
        for c = 1:backEndHandle.nClusters
            switch backEndHandle.statusRaw(c,1)
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
        display = cellfun(@(x) ~strcmp(x(1:3),'Con'),status,'uni',false);
        contam = num2cell(backEndHandle.contaminationValues);
        spcount = num2cell(backEndHandle.spikeCounts);
        sprate = num2cell(backEndHandle.spikeCounts * 20000 / backEndHandle.nSamples);
        comment = backEndHandle.comment;
        % Define the data
        % ID - Color - Display - Status - Contam
        d = [IDs, colors, display, status, spcount, sprate, contam, comment];
    end
    
    function selectorCallback(source,callbackdata)
        if source == allButton
            fakeCallbackData.NewData = true;
        else
            fakeCallbackData.NewData = false;
        end
        fakeCallbackData.Indices = [0,3];
        for c = 1:(backEndHandle.nClusters-1)
            clustMgmt.Data{c,3} = fakeCallbackData.NewData;
            fakeCallbackData.Indices(1) = c;
            tableEditCallback(source,fakeCallbackData,false);
        end
        c = backEndHandle.nClusters;
        clustMgmt.Data{c,3} = fakeCallbackData.NewData;
        fakeCallbackData.Indices(1) = c;
        tableEditCallback(source,fakeCallbackData,true);
    end
    
    function tableEditCallback(source,callbackdata,varargin)
        if ~(callbackdata.Indices(2) == 3)
            throw(MException('','ClusterEditGUI:tableEditCallback - Call from not a tickbox'));
        end
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
    
    function view3DCallback(source,callbackdata)
        if source == buttonDefaultView
            plot3D.View = [45 15];
        elseif source == buttonXY
            plot3D.View = [0 90];
        elseif source == buttonXZ
            plot3D.View = [0 0];
        elseif source == buttonYZ
            plot3D.View = [90 0];
        end
    end
    
    function autoScalePCPlots(varargin)
        listObj = [plot3D;PC45Plot];
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
    
    function s = bool2onoff(b)
        if b == 0
            s = 'off';
        else if b == 1
                s = 'on';
            else
                throw(MException('','bool2onoff: not 0 or 1'));
            end
        end
    end
    
    function s = bool2tf(b)
        if b == 0
            s = 'false';
        else if b == 1
                s = 'true';
            else
                throw(MException('','bool2tf: not 0 or 1'));
            end
        end
    end
    
    function [k,c] = nClustSelected()
        k = sum(cell2mat(clustMgmt.Data(:,3)));
        c = find(cell2mat(clustMgmt.Data(:,3)));
    end
    
    function customizeTools(figureHandle)
        allTools = findall(figureHandle);
        %%
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
    end
    
    frontEndHandle.Visible = 'on';
end