%function [backEndHandle,frontEndHandle] = ClusterEditGUI(datasetFolder)
function ClusterEditGUI
    %%
    mainFigure = figure( 'Name', 'Cluster Editor 0.0', ...
    'MenuBar', 'none', ...
    'Toolbar', 'figure', ...
    'NumberTitle', 'off',...
    'Visible', 'off');
    
    spacerWidth = 6;
    
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
    mainLayout.Sizes = [-2 -3];
    
    % graph box children
    % Graph 3D
    supportPanel3D = uipanel('Parent',graphLayout);
    plot3D = axes(...
        'Parent',supportPanel3D,...
        'ClippingStyle','rectangle');
    surf(plot3D,peaks);
    plot3D.Title.String = 'Principal Components 1-2-3';
    plot3D.XLabel.String = 'PC 1';
    plot3D.YLabel.String = 'PC 2';
    plot3D.ZLabel.String = 'PC 3';
    
    bottomGraphs = uiextras.HBox(...
        'Parent',graphLayout,...
        'Spacing',spacerWidth);
    graphLayout.Sizes = [-3 -2];
    
    % Two 2D plots at bottom right
    % Panels
    ACFPlot = axes(...
        'Parent',bottomGraphs);
    a = -5:0.01:5;
    plot(ACFPlot,a,a.^2,'k-+');
    ACFPlot.Title.String = 'ACF of selected clusters';
    ACFPlot.XLabel.String = 'Delay';
    ACFPlot.YLabel.String = 'Value';
    
    % 4th col bottom
    PC45Plot = axes(...
        'Parent',bottomGraphs);
    a = -5:0.01:5;
    plot(PC45Plot,a,sinc(a),'r--');
    PC45Plot.Title.String = 'Principal Components 4-5';
    PC45Plot.XLabel.String = 'PC 4';
    PC45Plot.YLabel.String = 'PC 5';
    
    bottomGraphs.Sizes = [-1 -1];
    
    datasetName = uicontrol(...
        'Parent',leftColumns,...
        'Style','text',...
        'fontsize',12,...
        'String','--The datasetName here--',...
        'Background','r');
    menuColumns = uiextras.HBox(...
        'Parent',leftColumns,...
        'Spacing',spacerWidth,...
        'Background','y');
    leftColumns.Sizes = [24 -1];
    
    % Menu Columns
    % Two hboxes
    leftMenu = uiextras.VBox(...
        'Parent',menuColumns,...
        'Spacing',spacerWidth,...
        'Background','g');
    rightMenu = uiextras.VBox(...
        'Parent',menuColumns,...
        'Spacing',spacerWidth,...
        'Background','b');
    menuColumns.Sizes = [-1 -1];
    
    % Left column menu - described by rows.
    loadRow = uiextras.HBox(...
        'Parent',leftMenu,...
        'Spacing',2*spacerWidth,...
        'Padding',spacerWidth);
    
    elNumberBox = uicontrol(...
        'Parent',loadRow,...
        'Style', 'edit',...
        'Max',1,'Min',1,...
        'String','El#',...
        'fontsize',11,...
        'callBack',@loadButtonCallback);
    
    loadButton = uicontrol(...
        'Parent',loadRow,...
        'Style', 'pushbutton',...
        'String','Load',...
        'fontsize',11,...
        'callback',@loadButtonCallback);
    
    fillerLoadRow = uicontrol(...
        'Parent',loadRow,...
        'Style','text',...
        'String','');
    
    clustNumberBox = uicontrol(...
        'Parent',loadRow,...
        'Style', 'edit',...
        'Max',1,'Min',1,...
        'String','Clust#',...
        'fontsize',11,...
        'callBack',@loadButtonCallback);
    
    loadClusterButton = uicontrol(...
        'Parent',loadRow,...
        'Style', 'pushbutton',...
        'String','Load',...
        'fontsize',11,...
        'callback',@loadClustButtonCallback);
    
    % finish loadRow
    loadRow.Sizes = [34, 50, -1,45, 60];
    
    % cluster setup table
    [columnname, columnformat, d] = getBackendDataClusterData();
    
    clustMgmt = uitable(...
        'Parent',leftMenu,...
        'Data', d,... 
        'ColumnName', columnname,...
        'ColumnFormat', columnformat,...
        'ColumnEditable', [false false true true],...
        'RowName',[],...
        'CellEditCallback',@testCellEditCallback);
    
    dummyButton = uicontrol(...
        'Parent',leftMenu,...
        'String','doStuff',...
        'Callback',@doStuff);
    
    
    % Finish left menu
    leftMenu.Sizes = [34 -1 20];
    
    function loadButtonCallback(source,callbackdata)
        fprintf('You clicked load! %s\n',elNumberBox.String);
        e = str2num(elNumberBox.String);
        if numel(e) ~= 1 || isnan(e)
            fprintf('Requested electrode not a number.\n');
            return;
        end
        % Check validity relative to dataset
        % call backend load
    end
    
    function loadClustButtonCallback(source,callbackdata)
        fprintf('You clicked load cluster! %s\n',clustNumberBox.String);
        c = str2num(elNumberBox.String);
        if numel(c) ~= 1 || isnan(c)
            fprintf('Requested cluster not a number.\n');
            return;
        end
        % Check validity relative to dataset
        % call backend load
    end
    
    function [columnname, columnformat, d] = getBackendDataClusterData()
        % HTML trick
        
    colorgen = @(color,text) ['<html><table border=0 width=400 bgcolor=',color,'><TR><TD>',text,'</TD></TR> </table></html>'];
    
        % Column names and column format
        columnname = {'Rate','Amount','Available','Fixed/Adj'};
        columnformat = {'numeric','char','logical',{'Fixed' 'Adjustable'}};
        % Define the data
        d ={6.125678  colorgen('#0000FF','A')  true   'Fixed';...
            6.75   colorgen('#00FF00','B')  false  'Adjustable';...
            7      colorgen('#FF0000','C')     false  'Fixed';};
    end
    
    function doStuff(source,callbackdata) % dummy callback table dynamic update
        clustMgmt.Data{1,1} = clustMgmt.Data{1,1} + 1;
    end
    
    function testCellEditCallback(source,callbackdata)
        source.Data{2,1} = source.Data{2,1} + 1;
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