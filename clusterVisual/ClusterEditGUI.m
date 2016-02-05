function [backEndHandle,frontEndHandle] = ClusterEditGUI(datasetFolder,varargin)
    % function ClusterEditGUI
    %% Instantiate - return handles
    narginchk(1,2);
    if nargin == 1
        backEndHandle = ClusterEditBackend(datasetFolder);
    else
        backEndHandle = varargin{1};
    end
    frontEndHandle = figure( 'Name', 'Cluster Editor 0.0', ...
        'MenuBar', 'none', ...
        'Toolbar', 'figure', ...
        'NumberTitle', 'off',...
        'Visible', 'off');
    
    %%
    mainFigure = frontEndHandle;
    
    spacerWidth = 10;
    displayPoints = 8000;
    
    clusterColors = zeros(0,3);
    
    % Main box and childre
    mainLayout = uiextras.HBoxFlex(...
        'Parent',mainFigure,...
        'Position',[0 0 1 1],...
        'Spacing',2);
    % Left VBox
    leftColumns = uiextras.VBox(...
        'Parent',mainLayout,...
        'Spacing',spacerWidth);
    
    graphLayout = uiextras.VBoxFlex(...
        'Parent',mainLayout,...
        'Spacing',spacerWidth);
    mainLayout.Sizes = [-1 -1];
    
    % graph box children
    % Graph 3D
    graph3DBox = uiextras.HBox(...
        'Parent',graphLayout,...
        'Spacing',spacerWidth);
    buttonColumn = uiextras.VBox(...
        'Parent',graph3DBox,...
        'Spacing',spacerWidth,...
        'Background','y');
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
    
    supportPanel3D = uipanel('Parent',graph3DBox);
    graph3DBox.Sizes = [50 -1];
    
    PC123Plots = {}; % handle to subplots handle defined in load callback
    plot3D = axes(...
        'Parent',supportPanel3D,...
        'ClippingStyle','rectangle',...
        'View',[45,15],...
        'XGrid','on','YGrid','on','ZGrid','on',...
        'DataAspectRatio',[1 1 1]);
    plot3D.Title.String = 'Principal Components 1-2-3';
    plot3D.XLabel.String = 'PC 1';
    plot3D.YLabel.String = 'PC 2';
    plot3D.ZLabel.String = 'PC 3';
    
    bottomGraphs = uiextras.HBox(...
        'Parent',graphLayout,...
        'Spacing',0*spacerWidth);
    graphLayout.Sizes = [-3 -2];
    
    % Three 2D plots at bottom right
    % Panels
    ratePlots = {}; % global handle to subplots handle defined in load callback
    ratePlot = axes(...
        'Parent',bottomGraphs,...
        'XGrid','on','YGrid','on');
    ratePlot.Title.String = 'Spike rates';
    ratePlot.XLabel.String = 'Time (sec)';
    ratePlot.YLabel.String = 'Spike rate (Hz)';
    ratePlot.XLim = [0, backEndHandle.nSamples / 20000 + 1];
    
    ACFPlots = {}; % global handle to subplots handle defined in load callback
    ACFPlot = axes(...
        'Parent',bottomGraphs,...
        'XGrid','on','YGrid','on');
    ACFPlot.Title.String = 'ACF';
    ACFPlot.XLabel.String = 'Time \Delta (msec)';
    ACFPlot.YLabel.String = 'Autocorr. (pair fraction/msec \Delta)';
    ACFPlot.XLim = [0, 101];
    
    % 4th col bottom
    PC45Plots = {}; % global handle to subplots handle defined in load callback
    PC45Plot = axes(...
        'Parent',bottomGraphs,...
        'XGrid','on','YGrid','on');
    PC45Plot.Title.String = 'Principal Components 4-5';
    PC45Plot.XLabel.String = 'PC 4';
    PC45Plot.YLabel.String = 'PC 5';
    
    bottomGraphs.Sizes = [-1 -1 -1];
    
    datasetName = uicontrol(...
        'Parent',leftColumns,...
        'Style','text',...
        'fontsize',12,...
        'String',datasetFolder);
    menu = uiextras.VBox(...
        'Parent',leftColumns,...
        'Spacing',0,...
        'Background','g');
    eistaBox = uiextras.HBox('Parent',leftColumns,...
        'Spacing',0,'background','b');
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
        'String','El#',...
        'fontsize',11,...
        'callback',@loadButtonCallback);
    
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
        'String','',...
        'Background',[0 0.5 0.5]);
    
    % finish loadRow
    loadRow.Sizes = [34, 50, 45, 60, -1];
    
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
    refreshButton = uicontrol(...
        'Parent',selectorRow,...
        'Style', 'pushbutton',...
        'String','Re-scatter',...
        'fontsize',11,...
        'callback',@refreshGraphics);
    fillerSelectorRow = uicontrol(...
        'Parent',selectorRow,...
        'Style','text',...
        'String','',...
        'Background',[0.5 0.5 0]);
    
    selectorRow.Sizes = [75 90 75 -1];
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
    
    eiPanel = uipanel('Parent',eistaBox);
    eiJPanel = []; eiHPanel = []; % Globalize the java panels
    eiDistAxes = axes('Parent',eiPanel,...
        'Visible','off');
    noEIString = uicontrol(...
        'Parent',eiPanel,...
        'Style','text',...
        'fontsize',12,...
        'String','No EI for this neuron',...
        'Visible','off');
    noEIString.Units = 'Norm';
    noEIString.Position = [0 0 1 1];
    
    staPanel = uipanel('Parent',eistaBox);
    staJPanel = []; staHPanel = []; % Globalize the java panels
    corrDistAxes = axes('Parent',staPanel,...
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
    function loadButtonCallback(source,callbackdata)
        statusBar.String = sprintf('Loading electrode %s...',elNumberBox.String);
        e = str2num(elNumberBox.String);
        if numel(e) ~= 1 || isnan(e)
            statusBar.String = sprintf('Requested electrode %s not a number.',elNumberBox.String);
            return;
        end
        status = backEndHandle.loadEl(e);
        if status ~= 0
            return;
        end
        refreshView();
        if source == loadButton || source == elNumberBox
            clustNumberBox.String = 'ID#';
        end
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
            elNumberBox.String = el;
            loadButtonCallback(source,callbackdata);
        end
    end
    
    function refreshView()
        makeColors();
        clustMgmt.Data = getClusterData(); pause(.1);
        refreshGraphics();
        refreshLowLeftPanels();
    end
    
    function makeColors()
        % make HSV colors
        %clusterColors.HSV = [randperm(backEndHandle.nClusters)' / backEndHandle.nClusters,...
        %    ones(backEndHandle.nClusters,1)*[0.8 0.7]];
        clusterColors.HSV = [(1:backEndHandle.nClusters)' / backEndHandle.nClusters,...
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
            % EI Panel
            if numel(backEndHandle.eiFile) > 0 && numel(backEndHandle.eisLoaded{c}) > 0
                eiDistAxes.Visible = 'off';
                noEIString.Visible = 'off';
                [eiJPanel,eiHPanel] = ...
                    javacomponent(edu.ucsc.neurobiology.vision.neuronviewer.PhysiologicalImagePanel(...
                    backEndHandle.eisLoaded{c},...
                    [],2,...
                    backEndHandle.electrodeMap,...
                    backEndHandle.elLoaded - 1));
                eiHPanel.Parent = eiPanel;
                eiHPanel.Units = 'norm';
                eiHPanel.Position = [0 0 1 1];
                isEI = true;
            else
                isEI = false;
                eiDistAxes.Visible = 'off';
                noEIString.Visible = 'on';
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
                corrDistAxes.Visible = 'off';
                noSTAString.Visible = 'on';
            end
            
            if ~isEI && ~isSTA
                leftColumns.Sizes(3) = 26;
            end
            
        else % Multiple neurons selected - display EI distance and xCorr distance
            leftColumns.Sizes(3) = -1;
            % EI Panel (Left)
            noEIString.Visible = 'off';
            
            eiMatrix = imagesc(backEndHandle.EIdistMatrix,...
                'AlphaData',~isnan(backEndHandle.EIdistMatrix),...
                'Parent',eiDistAxes);
            title(eiDistAxes,'Pairwise EI distance');
            eiDistAxes.DataAspectRatio = [1 1 1];
            eiDistAxes.XAxisLocation = 'top';
            eiDistAxes.XTick = 1:backEndHandle.nClusters;
            eiDistAxes.YTick = 1:backEndHandle.nClusters;
            eiDistAxes.XTickLabel = num2cell(backEndHandle.displayIDs);
            eiDistAxes.YTickLabel = num2cell(backEndHandle.displayIDs);
            eiDistAxes.XTickLabelRotation = 35;    
            eiDistAxes.YTickLabelRotation = 35;
            eiDistAxes.TickLength = [0,0];
            colorbar(eiDistAxes,'eastoutside');
            eiDistAxes.Position = [0.1    0    0.68    0.9];
            eiDistAxes.CLim = [0,1];
            eiDistAxes.Visible = 'on';
            
            nc = backEndHandle.nClusters;
            hold(eiDistAxes,'on');
            plot(eiDistAxes,repmat([0.5,nc + 0.5],nc-1,1)',((0.5+(1:(nc-1)))'*[1,1])','k-','linewidth',1);
            plot(eiDistAxes,((0.5+(1:(nc-1)))'*[1,1])',repmat([0.5,nc + 0.5],nc-1,1)','k-','linewidth',1);
            hold(eiDistAxes,'off');
            
            % STA Panel (Right)
            noSTAString.Visible = 'off';
            
            corrMatrix = imagesc(backEndHandle.spikeTrainCorr,...
                'AlphaData',~isnan(backEndHandle.spikeTrainCorr),...
                'Parent',corrDistAxes);
            title(corrDistAxes,'Spike train correlation');
            colormap(corrDistAxes,flipud(colormap(eiDistAxes)));
            corrDistAxes.DataAspectRatio = [1 1 1];
            corrDistAxes.XAxisLocation = 'top';
            corrDistAxes.XTick = 1:backEndHandle.nClusters;
            corrDistAxes.YTick = 1:backEndHandle.nClusters;
            corrDistAxes.XTickLabel = num2cell(backEndHandle.displayIDs);
            corrDistAxes.YTickLabel = num2cell(backEndHandle.displayIDs);
            corrDistAxes.XTickLabelRotation = 35;    
            corrDistAxes.YTickLabelRotation = 35;
            corrDistAxes.TickLength = [0,0];
            colorbar(corrDistAxes,'eastoutside');
            corrDistAxes.Position = [0.1    0    0.68    0.9];
            corrDistAxes.CLim = [0,1];
            corrDistAxes.Visible = 'on';
            
            nc = backEndHandle.nClusters;
            hold(corrDistAxes,'on');
            plot(corrDistAxes,repmat([0.5,nc + 0.5],nc-1,1)',((0.5+(1:(nc-1)))'*[1,1])','k-','linewidth',1);
            plot(corrDistAxes,((0.5+(1:(nc-1)))'*[1,1])',repmat([0.5,nc + 0.5],nc-1,1)','k-','linewidth',1);
            hold(corrDistAxes,'off');
        end
        % eiPanel.Children
        % staPanel.Children
    end
    
    % Displays the info table in the left columns
    % operates by side effects
    % To split:
    % Column headers and format when instantiating the table
    % table data when refreshing after load.
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
        status = backEndHandle.status;
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
    
    
    
    
    
    % Instantiate backend
    % Backend deals with all argument checking, partial arguments, etc...
    % And data management
    
    % Instantiate figure
    % Instatiate axes/subplots
    % Plot PC 1-2-3
    % Plot PC 4-5
    % Plot ACF
    %
    
    % Instantiate buttons
    % Each button:
    % Type
    % Position
    % Specifics
    % Callback
    
    % For plot PC 1-2-3
    % x view num boxes
    % y view num boxes
    % z view num boxes
    % "Refresh view"
    
    % For plot PC 4-5
    % x view boxes
    % y view boxes
    
    mainFigure.Visible = 'on';
    'done'
end