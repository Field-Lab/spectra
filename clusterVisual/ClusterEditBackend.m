classdef ClusterEditBackend < handle
    %CLUSTEREDITBACKEND Handles persistently the backend of the cluster visualizer
    
    properties (GetAccess = public, SetAccess = private)
        typeIsSpectra
        
        analysisPath
        
        prj % struct - {bool exists ; String path ; matfile or javafile}
        neurons % struct - {bool exists ; String path ; matfile or javafile}
        
        config
        
        % Java related stuff for ei and sta display
        eiFileRaw
        eiFile
        staFileRaw
        staFile
        globalsFile
        
        arrayID
        electrodeMap
        
        savedState
    end
    
    properties(GetAccess = public,SetAccess = protected)
        % General dataset information
        nNeurons % number of neurons
        nElectrodes % number of electrodes
        nSamples % number of samples
        
        neuronEls % electrode of each neuron
        neuronClusters % cluster number of each cluster
        neuronIDs % neuron ID of each cluster
        neuronStatuses % cleaning status report for each cluster
        classification % classification string for each cluster
        
        elLoaded % electrode currently loaded - MATLAB numbering
        
        
        % Front end friendly data
        isDataReady % Marker if loading is complete. Is not currently used
        
        % Information about the clusters on loaded electrode
        nClusters % Number of clusters
        displayIDs % IDs
        contaminationValues % Contamination values
        spikeCounts % Spike train counts
        statusRaw % Raw clean status markers
        comment % string parsed from classification file
        
        spikeTrains % spike trains of clusters
        prjTrains % Projections of the spike train
        eisLoaded % (available) EIs
        EIdistMatrix % computed EI distance matrix
        stasLoaded % (available) STAs
    end
    
    properties % public setting access
        statusBarHandle = [] % Handle to the GUI status bar - set from the GUI
        % To keep working without a GUI, a statusBarHandle.String = '' initializer is in the
        % constructor
    end
    
    methods
        % Constructor
        % Builds backend object defined in this file for the cluster visualizer
        % Input
        %       analysisPath: path to analysis folder
        %       OPTIONAL (none or all three)
        %       pathToPrj: full path of .prj.mat file to use
        %       pathToNeurRaw: full path to .neurons.mat to use
        %       pathToModel: full path to .model.mat to use (unused so far)
        %
        % One or more of the optional arguments can be set to empty string ''
        % in which case constructor looks for a file with correct extension
        % in the analysis path.
        function obj = ClusterEditBackend(analysisPath)
            obj.analysisPath = analysisPath;
            obj.config = mVisionConfig();
            
            % Argument check
            validateattributes(analysisPath,{'char'},{},'','Analysis Path',1);
            if ~ (exist(analysisPath,'file') == 7)
                throw(MException('','ClusterEditBackend:ClusterEditBackend - Analysis folder does not exist'));
            end
            
            % Define the type
            files = dir([analysisPath,filesep,'*.prj.mat']);
            if numel(files) == 1
                obj.typeIsSpectra = 1;
                fprintf('Starting in Spectra mode.\n');
            else
                obj.typeIsSpectra = 0;
                fprintf('Starting in Vision mode. Edition features disabled.\n')
            end
            
            % Type specific loading
            if obj.typeIsSpectra
                % Projections file
                files = dir([analysisPath,filesep,'*.prj.mat']);
                if numel(files) == 1
                    fprintf('Loading projections file...\n');
                    obj.prj.exist = true;
                    obj.prj.path = [analysisPath,filesep,files(1).name];
                else if numel(files) == 0
                        throw(MException('','ClusterEditBackend:ClusterEditBackend - No projections file found'));
                    else
                        throw(MException('','ClusterEditBackend:ClusterEditBackend - Multiple projections files found\nPlease use extended constructor call\nneuronViewer(analysisPath,projectionsFile,neuronsRawFile,modelFile)'));
                    end
                end
                
                % Neurons file
                files = dir([analysisPath,filesep,'*-edited.neurons.mat']);
                if numel(files) == 0 % No edited - start from raw
                    files = dir([analysisPath,filesep,'*.neurons.mat']);
                    fprintf('Loading neurons from raw neurons file.\n')
                else
                    fprintf('Loading clusters from consolidated neurons file.\n');
                end
                if numel(files) == 1
                    obj.neurons.exist = true;
                    obj.neurons.path = [analysisPath,filesep,files(1).name];
                else if numel(files) == 0
                        throw(MException('','ClusterEditBackend:ClusterEditBackend - No neurons.mat file found'));
                    else
                        throw(MException('','ClusterEditBackend:ClusterEditBackend - Multiple neurons.mat files found\nPlease use extended constructor call\nneuronViewer(analysisPath,projectionsFile,neuronsRawFile,modelFile)'));
                    end
                end
                
                % All files now referenced and checked
                % Partial loading for neuronSpikeTimes and elSpikeTimes
                obj.prj.matfile = matfile(obj.prj.path);
                obj.neurons.matfile = matfile(obj.neurons.path);
                
                load(obj.neurons.path,'neuronEls','neuronClusters','nSamples'); % Loads neuronClusters, neuronEls
                obj.nSamples = nSamples;
                obj.neuronEls = neuronEls;
                obj.neuronClusters = neuronClusters;
                obj.neuronIDs = ClusterEditBackend.getIDs(obj.neuronEls, obj.neuronClusters);
                
                obj.nElectrodes = size(obj.prj.matfile.spikeTimes,1);
                obj.nNeurons = size(obj.neuronEls,1);
                
                % Clean pattern %
                files = dir([analysisPath,filesep,'*.clean.mat']);
                obj.neuronStatuses = zeros(obj.nNeurons,2);
                % check against elevated statuses
                if numel(files) == 1
                    fprintf('Parsing automatic action stack from .clean.mat file.\n');
                    % -2 - discard ; -1 - Unknown ; 0 - keep ; 1 - removed ; 2 - merged (2nd col: master); 3 - dup (2nd col: master)
                    % Parse statuses
                    % Avoid interference with elevated IDs
                    load(obj.neurons.path,'elevatedStatus');
                                        if (~exist('elevatedStatus', 'var')) 
                        elevatedStatus = false(obj.nNeurons,1);
                    end
                    listElevatedRows = find(elevatedStatus);
                    
                    % Fast positional neuron access. fastRows(ID) = row #.
                    fastRows = zeros(max(obj.neuronIDs) + 15,1);
                    % + 15 safety in case manual actions add IDs
                    [~,~,positions] = intersect(1:max(obj.neuronIDs),obj.neuronIDs);
                    fastRows(obj.neuronIDs) = positions;
                    
                    cleanPatternPath = [analysisPath,filesep,files(1).name];
                    load(cleanPatternPath); % loads autoActions
                    
                    % Obtain status from dry-running all auto cleaning actions
                    % Apply all automatic actions
                    for i = 1:size(autoActions,1)
                        action = autoActions{i,1};
                        params = autoActions{i,2};
                        switch action
                            case EditAction.AUTO_RM
                                r = fastRows(params{1}); r(r == 0) = []; % find the rows
                                r = setdiff(r,listElevatedRows); % remove the elevated rows
                                obj.neuronStatuses(r,1) = 1; % Set status to 1
                            case EditAction.AUTO_MERGE
                                rMaster = fastRows(params{1}); % find the rows
                                r = fastRows(params{2}); r(r == 0) = []; % remove elevated rows
                                r = setdiff(r,listElevatedRows);
                                % depends on whether master ID is elevated
                                if ~(rMaster == 0) && ~elevatedStatus(rMaster)
                                    obj.neuronStatuses(r,1) = 2;
                                    obj.neuronStatuses(r,2) = params{1}; % Master ID
                                else % master is elevated, need to find another master neuron
                                    if numel(r) > 0
                                        nSpikes = cellfun(@(x) numel(x),obj.neurons.matfile.neuronSpikeTimes(r),'uni',true);
                                        [~,b] = max(nSpikes);
                                        rMaster = r(b);
                                        r(b) = [];
                                        obj.neuronStatuses(r,1) = 2;
                                        obj.neuronStatuses(r,2) = obj.neuronIDs(rMaster);
                                    end
                                end
                            case EditAction.AUTO_RM_DUP
                                r = fastRows(params{2});  r(r == 0) = []; % find the rows
                                r = setdiff(r,listElevatedRows); % remove elevated rows
                                obj.neuronStatuses(r,1) = 3;
                                obj.neuronStatuses(r,2) = params{1};
                            otherwise
                                throw(MException('','ClusterEditBackend:ClusterEditBackend - Unhandled EditAction in switch statement.'));
                        end
                    end
                    % shorten duplicate chains
                    dupRows = obj.neuronStatuses(:,1) == 3;
                    tmp = ClusterEditBackend.shortenDuplicatesPath([obj.neuronStatuses(dupRows,2),obj.neuronIDs(dupRows)]);
                    obj.neuronStatuses(dupRows,2) = tmp(:,1);
                else
                    fprintf('No .clean.mat file (or more than one). Skipping automatic action stack parsing.\n');
                end
            else % Type is vision
                % Projections file
                files = dir([analysisPath,filesep,'*.prj']);
                if numel(files) == 1
                    fprintf('Loading projections file...\n');
                    obj.prj.exist = true;
                    obj.prj.path = [analysisPath,filesep,files(1).name];
                else if numel(files) == 0
                        throw(MException('','ClusterEditBackend:ClusterEditBackend - No projections file found'));
                    else
                        throw(MException('','ClusterEditBackend:ClusterEditBackend - Multiple projections files found\nPlease use extended constructor call\nneuronViewer(analysisPath,projectionsFile,neuronsRawFile,modelFile)'));
                    end
                end
                
                % Neurons file
                files = dir([analysisPath,filesep,'*.neurons-raw']);
                if numel(files) == 1
                    fprintf('Loading neurons file...\n');
                    obj.neurons.exist = true;
                    obj.neurons.path = [analysisPath,filesep,files(1).name];
                else if numel(files) == 0
                        throw(MException('','ClusterEditBackend:ClusterEditBackend - No .neurons-raw file found'));
                    else
                        throw(MException('','ClusterEditBackend:ClusterEditBackend - Multiple .neurons-raw files found\nPlease use extended constructor call\nneuronViewer(analysisPath,projectionsFile,neuronsRawFile,modelFile)'));
                    end
                end
                
                % All files now referenced and checked
                % Partial loading for neuronSpikeTimes and elSpikeTimes
                obj.prj.javaFile = edu.ucsc.neurobiology.vision.io.ProjectionsFile(obj.prj.path);
                obj.neurons.javaFile = edu.ucsc.neurobiology.vision.io.NeuronFile(obj.neurons.path);
                
                % Generate information
                
                hdr = obj.neurons.javaFile.getHeader();
                obj.nSamples = hdr.nSamples;
                obj.neuronIDs = double(obj.neurons.javaFile.getIDList());
                [obj.neuronEls,obj.neuronClusters] = ClusterEditBackend.getElClust(obj.neuronIDs);
                obj.nNeurons = size(obj.neuronEls,1);
                % Correct way - but header contains garbage arrayID
                % elMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(hdr.arrayID);
                % obj.nElectrodes = elMap.getNumberOfElectrodes();
                obj.nElectrodes = max(obj.neuronEls); % Bypass
                % no clean pattern in vision mode, unless we can find a final neurons file
                obj.neuronStatuses = zeros(obj.nNeurons,2);
                files = dir([analysisPath,filesep,'*.neurons']);
                if numel(files) == 1
                    fprintf('Parsing neuron statuses from post-cleaning neurons file.\n');
                    neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile([analysisPath,filesep,files(1).name]);
                    IDsKept = double(neuronFile.getIDList());
                    neuronFile.close();
                    [~,i,~] = intersect(obj.neuronIDs,IDsKept);
                    obj.neuronStatuses(:,1) = -2;
                    obj.neuronStatuses(i,1) = 0;
                else
                    obj.neuronStatuses(:,1) = -1;
                    fprintf('No post-cleaning neurons file (or more than one), neuron statuses are unknown.\n');
                end
            end % type switch (end vision type)
            
            % Initialize EI
            obj.initEI(analysisPath);
            
            % Initialize STA
            obj.initSTA(analysisPath);
            
            % Try to find a classification file
            obj.classification = cell(obj.nNeurons,1);
            files = dir([analysisPath,filesep,'*.classification.txt']);
            if numel(files) == 1
                fprintf('Loading labels from .classification.txt file.\n');
                fid = fopen([analysisPath,filesep,files(1).name]);
                classesRaw = textscan(fid, '%u All/%s', 'delimiter', '\n');
                [~,pos,pos2] = intersect(obj.neuronIDs,classesRaw{1});
                obj.classification(pos) = classesRaw{2}(pos2);
                fclose(fid);
            else
                fprintf('No .classification.txt file found (or multiple). Labels unavailable.\n');
            end
            
            % electrode map
            if numel(obj.eiFile) > 0
                obj.arrayID = obj.eiFile.arrayID;
                obj.electrodeMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(obj.arrayID);
            end
            
            % Handling case where no GUI caller exists - user mode
            % Code will pass through calls to the (fake) GUI status bar no questions asked.
            obj.statusBarHandle.String = '';
            
            fprintf('All data loaded - ClusterEdit backend is ready.\n');
        end
        
        % function initEI
        %   Relocation of a piece of the constructor
        %   Setting up ei and ei-raw access
        function initEI(obj, analysisPath)
            % Initialize EI - RAWS
            files = dir([analysisPath,filesep,'*.ei-raw']);
            if numel(files) > 0
                fprintf('Loading raw EI file...\n');
                eiPath = [analysisPath,filesep,files(1).name];
                obj.eiFileRaw = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(eiPath);
            else
                obj.eiFileRaw = []; % Tag for existence of ei File
                fprintf('No EI (RAW) file found, skipping\n');
            end
            
            % Initialize EI - FINAL
            files = dir([analysisPath,filesep,'*.ei']);
            if numel(files) > 0
                fprintf('Loading final EI file...\n');
                eiPath = [analysisPath,filesep,files(1).name];
                obj.eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(eiPath);
            else
                obj.eiFile = []; % Tag for existence of ei File
                fprintf('No EI (FINAL) file found, skipping\n');
            end
        end
        
        
        % function initSTA
        %   Relocation of a piece of the constructor
        %   Setting up sta and sta-raw access
        function initSTA(obj, analysisPath)
            files = dir([analysisPath,filesep,'*.sta-raw']);
            if numel(files) > 0
                fprintf('Loading raw STA file...\n');
                staPath = [analysisPath,filesep,files(1).name];
                obj.staFileRaw = edu.ucsc.neurobiology.vision.io.STAFile(staPath);
            else
                obj.staFileRaw = []; % Tag for existence of ei File
                fprintf('No STA (RAW) file found, skipping\n');
            end
            
            files = dir([analysisPath,filesep,'*.sta']);
            if numel(files) > 0
                fprintf('Loading final STA file...\n');
                staPath = [analysisPath,filesep,files(1).name];
                obj.staFile = edu.ucsc.neurobiology.vision.io.STAFile(staPath);
            else
                obj.staFile = []; % Tag for existence of ei File
                fprintf('No STA (FINAL) file found, skipping\n');
            end
            
            % Grab the globals file
            files = dir([analysisPath,filesep,'*.globals']);
            if numel(files) > 0
                globalsPath = [analysisPath,filesep,files(1).name];
                obj.globalsFile = edu.ucsc.neurobiology.vision.io.chunk.GlobalsFile(globalsPath);
            else
                obj.globalsFile = [];
                obj.staFile = [];
                obj.staFileRaw = [];
                fprintf('ClusterEditBackend:ClusterEditBackend - No Globals file found, STAs Unavailable\n');
            end
        end
        
        % function loadEl
        %   loads the requested electrode and prepares all front-end friendly data
        %   Input:
        %       el: MATLAB numbered electrode number
        %       optional varargin{1}: put to 1 to override behavior that prevent reloading data
        %           already loaded.
        %   Output:
        %       returnStatus: 1 if load failed (same as already loaded, invalid)
        %                     0 if load succeded
        function returnStatus = loadEl(obj,el,varargin)
            
            if numel(varargin) == 0 || ~varargin{1} % Already loaded skip override
                if el == obj.elLoaded % Already loaded - save the work
                    obj.statusBarHandle.String = sprintf('Electrode %u already loaded.',el-1);
                    returnStatus = 1;
                    return;
                end
            end
            if (el <= 1) || (el > obj.nElectrodes) % Invalid number
                obj.statusBarHandle.String = sprintf('Electrode number %u invalid, nothing done.',el-1);
                returnStatus = 1;
                return;
            end
            
            % Clear the restore point
            obj.saveLocalState('clear');
            
            % Load the electrode
            obj.elLoaded = el;
            if obj.typeIsSpectra % Matfile prj loading in spectra mode
                eval(sprintf('load(''%s'',''projSpikes%u'');',obj.prj.path,el));
                eval(sprintf('prjLoaded = projSpikes%u;',el));
            else % Java loading in vision mode
                projToMatlabWrap = obj.prj.javaFile.readProjections(el-1);
                prjLoaded = double(projToMatlabWrap.data)';
            end
            % position of relevant neurons in the global list
            neuronIndices = find(obj.neuronEls == el);
            
            % Electrode data - initializers
            obj.displayIDs = obj.neuronIDs(neuronIndices);
            obj.nClusters = numel(obj.displayIDs);
            obj.spikeTrains = cell(obj.nClusters,1);
            obj.prjTrains = cell(obj.nClusters,1);
            obj.contaminationValues = zeros(obj.nClusters,1);
            obj.spikeCounts = zeros(obj.nClusters,1);
            obj.statusRaw = obj.neuronStatuses(neuronIndices,:);
            obj.comment = obj.classification(neuronIndices);
            
            if obj.typeIsSpectra
                % Get all spike times on electrode from prj.mat file
                % Reconstruct projection trains by cluster
                elSpikeTimes = obj.prj.matfile.spikeTimes(obj.elLoaded,1);
                elSpikeTimes = elSpikeTimes{1};
                
                if numel(neuronIndices) > 0
                    % Load spike trains
                    obj.spikeTrains = obj.neurons.matfile.neuronSpikeTimes(neuronIndices,1);
                else
                    obj.spikeTrains = cell(0,1);
                end
            else
                % Get all spike times on electrode from the prj file
                % projectionsToMatlab returned object
                % Reconstruct projections trains
                elSpikeTimes = double(projToMatlabWrap.times);
                obj.spikeTrains = cell(numel(obj.displayIDs),1);
                for c = 1:obj.nClusters
                    obj.spikeTrains{c} = double(obj.neurons.javaFile.getSpikeTimes(obj.displayIDs(c)));
                end
            end
            
            % Compute front-end ready data
            obj.spikeCounts = cellfun(@numel, obj.spikeTrains,'uni',true);
            
            for c = 1:obj.nClusters
                [~,~,indices] = intersect(...
                    obj.spikeTrains{c},...
                    elSpikeTimes);
                obj.prjTrains{c} = prjLoaded(indices,:);
                
                obj.contaminationValues(c) = ...
                    edu.ucsc.neurobiology.vision.anf.NeuronCleaning.getContam(obj.spikeTrains{c},int32(obj.nSamples));
            end
            
            % Load EIs
            if numel(obj.eiFile) > 0
                obj.loadEI();
                obj.loadEIDistances();
            else
                obj.EIdistMatrix = nan(obj.nClusters);
            end
            % Load STAs
            if numel(obj.staFile) > 0
                obj.loadSTA();
            end
            
            % Finish
            obj.isDataReady = true;
            returnStatus = 0;
        end
        
        % function reload same data
        %   reloads the electrode already selected from files and global dataset variables into
        %   local loaded electrode variables.
        function reloadSameData(obj)
            obj.loadEl(obj.elLoaded,true);
        end
        
        
        % function loadEI
        %   loads the EIs from the EI file for the IDs in obj.displayIDs
        function loadEI(obj)
            obj.eisLoaded = cell(obj.nClusters,1);
            for c = 1:obj.nClusters
                try
                    obj.eisLoaded{c} = obj.eiFile.getImage(obj.displayIDs(c));
                    if numel(obj.eisLoaded{c}) == 0
                        obj.eisLoaded{c} = obj.eiFileRaw.getImage(obj.displayIDs(c));
                    end
                catch
                    obj.eisLoaded{c} = [];
                end
            end
        end
        
        % function loadEIDistances
        %   Computes the EI distance matrix for the EIs loaded in obj.eisLoaded
        %   Metric identical to duplicate removal (realignement per electrode, L2 norm,...)
        function loadEIDistances(obj)
            cleanConfig = obj.config.getCleanConfig();
            eiSize = obj.eiFile.nlPoints + obj.eiFile.nrPoints + 1;
            
            % Recenter global minima window if ei left and right are not identical to the configuration
            cleanConfig.globMinWin = cleanConfig.globMinWin - cleanConfig.EILP + obj.eiFile.nlPoints;
            
            % Utils for realignment
            minWindowSize = cleanConfig.globMinWin(2) - cleanConfig.globMinWin(1);
            samplingVal = cleanConfig.globMinWin(1):cleanConfig.resamplePitch:cleanConfig.globMinWin(2);
            basicValues = 1:(eiSize - minWindowSize);
            
            % Electrode neighborhoods
            [adjacent2,~] = catchAdjWJava( obj.eiFile, 2);
            
            trace = zeros(obj.nClusters, (eiSize - minWindowSize) * numel(adjacent2{obj.elLoaded}));
            
            % EI upsampling handling
            traceTemp = zeros(eiSize - minWindowSize, numel(adjacent2{obj.elLoaded}));
            hasEI = cellfun(@(x) numel(x) > 0,obj.eisLoaded);
            for n = 1:obj.nClusters
                if ~hasEI(n)
                    continue
                end
                % Realignment procedure
                for neighborElInd = 1:numel(adjacent2{obj.elLoaded});
                    neighborEl = adjacent2{obj.elLoaded}(neighborElInd);
                    
                    interpolant = griddedInterpolant(1:eiSize,squeeze(obj.eisLoaded{n}(1,neighborEl,:)),'spline');
                    [~, offset] = min(interpolant(samplingVal));
                    traceTemp(:,neighborElInd) = interpolant(basicValues + (offset - 1) * cleanConfig.resamplePitch);
                end
                trace(n,:) = traceTemp(:)';
            end % n
            % L2 normalization
            trace(hasEI,:) = bsxfun(@rdivide,trace(hasEI,:),sqrt(sum(trace(hasEI,:).^2,2)));
            
            % Distance computation - nans for unavailable EIs
            [~,dists] = returnL2MergeClasses(trace(hasEI,:), 0);
            obj.EIdistMatrix = nan(obj.nClusters);
            obj.EIdistMatrix(hasEI,hasEI) = dists;
        end
        
        % Function loadSTA
        %   load the STAs in the STA file for the IDs in obj.displayIDs
        function loadSTA(obj)
            obj.stasLoaded = cell(obj.nClusters,1);
            for c = 1:obj.nClusters
                try
                    obj.stasLoaded{c} = obj.staFile.getSTA(obj.displayIDs(c));
                    % Try raw load
                    if numel(obj.stasLoaded{c}) == 0
                        obj.stasLoaded{c} = obj.staFileRaw.getSTA(obj.displayIDs(c));
                    end
                catch
                    obj.stasLoaded{c} = [];
                end
            end
        end
        
        % function checkID
        %   returns the electrode number for an ID of the dataset
        %   Input:
        %       ID: cluster ID
        %   Output:
        %       el: electrode number for this ID if it exists, -1 otherwise
        function el = checkID(obj,ID)
            rowNum = find(obj.neuronIDs == ID);
            if numel(rowNum) ~= 1
                obj.statusBarHandle.String = sprintf('Neuron ID invalid, nothing done.');
                el = -1;
                return;
            end
            [el,~] = obj.getElClust(ID);
        end
        
        % function localApplyAction
        %   applies an edit action to the data stored in the cache of this object,
        %   locally, on the selected electrode.
        %   This is previewing feature for edition actions.
        %
        %   Inputs:
        %       action: EditAction object
        %       parameters: cell array describing action parameters
        %       optional: data for fast reapply of shrink/recluster
        %
        %   Returns:
        %       returnStat: flag at 0 if OK, 1 if error
        %       data: data about the operation result, for pass to editHandler through the GUI
        %             error message to report in case of failure
        function [returnStat,data] = localApplyAction(obj,action,params,varargin)
            validateattributes(action,{'EditAction'},{});
            [v,m] = action.checkParameters(params);
            if v == 0
                throw(MException('',['ClusterEditBackend:localApplyAction - Invalid action parameters: \n',m]));
            end
            % Local copy of the previously stored previous state
            % Needs to be put back
            % obj.savedState = localStateCopy;
            % In case of any nonterminal return of this function
            localStateCopy = obj.savedState;
            obj.saveLocalState('save');
            
            switch action
                case EditAction.KEEP
                    [~,selRowsIdx,~] = intersect(obj.displayIDs,params{1});
                    if numel(selRowsIdx) < numel(params{1})
                        data = 'Some requested IDs are invalid. Aborting edition.';
                        returnStat = 1;
                        obj.savedState = localStateCopy;
                        return;
                    end
                    obj.statusRaw(selRowsIdx,:) = 0;
                    obj.comment(selRowsIdx) = cellfun(@(s) [s,' - Elevated'],obj.comment(selRowsIdx),'uni',false);
                    data = {};
                case EditAction.DELETE
                    [~,selRowsIdx,~] = intersect(obj.displayIDs,params{1});
                    if numel(selRowsIdx) < numel(params{1})
                        data = 'Some requested IDs are invalid. Aborting edition.';
                        returnStat = 1;
                        obj.savedState = localStateCopy;
                        return;
                    end
                    obj.statusRaw(selRowsIdx,:) = -2;
                    obj.comment(selRowsIdx) = cellfun(@(s) [s,' - Deleted'],obj.comment(selRowsIdx),'uni',false);
                    data = {};
                case EditAction.MERGE
                    [~,selRowsIdx,~] = intersect(obj.displayIDs,params{1});
                    if numel(selRowsIdx) < numel(params{1})
                        data = 'Some requested IDs are invalid. Aborting edition.';
                        returnStat = 1;
                        obj.savedState = localStateCopy;
                        return;
                    end
                    rMaster = find(obj.spikeCounts == max(obj.spikeCounts(selRowsIdx)),1); % row of final neuron in local space
                    saveCom = obj.comment{rMaster};
                    mergedTrain = sort(horzcat(obj.spikeTrains{selRowsIdx}));
                    % Use a reapply mode recluster action to process the merge.
                    obj.localApplyAction(EditAction.RECLUSTER,[params{1};{1};{''}],...
                        {obj.displayIDs(rMaster), {mergedTrain}});
                    % rMaster is now invalid - RECLUSTER put row at the end.
                    [~,rMaster,~] = intersect(obj.displayIDs,params{1});
                    obj.comment{rMaster} = [saveCom,' - Merged']; 
                    
                    contam = edu.ucsc.neurobiology.vision.anf.NeuronCleaning.getContam(mergedTrain,int32(obj.nSamples));
                    data = {numel(mergedTrain),... % #Spikes
                        numel(mergedTrain) * 20000./obj.nSamples,... % Spike rate
                        contam,...
                        mergedTrain }; % contamination value
                    
                case EditAction.RECLUSTER
                    % Switch cases whether applying stored action or actually computing
                    if numel(varargin) == 0
                        % First check that the config tag is valid
                        global GLOBAL_CONFIG
                        try
                            GLOBAL_CONFIG = mVisionConfig(params{3});
                        catch error
                            data = 'The configuration tag requested is invalid';
                            returnStat = 1;
                            obj.savedState = localStateCopy;
                            return;
                        end
                    end
                    % Second check ID validity
                    [~,selRowsIdx,~] = intersect(obj.displayIDs,params{1});
                    if numel(selRowsIdx) < numel(params{1})
                        data = 'Some requested IDs are invalid. Aborting edition.';
                        returnStat = 1;
                        obj.savedState = localStateCopy;
                        return;
                    end
                    
                    % Merged pool of spikes to recluster
                    allPrjs = vertcat(obj.prjTrains{selRowsIdx});
                    allTimes = horzcat(obj.spikeTrains{selRowsIdx});
                    [allTimes,sortTime] = sort(allTimes,'ascend');
                    allPrjs = allPrjs(sortTime,:);
                    
                    if numel(varargin) == 0 % Computing a newly requested action
                        try
                            [clusterIndexes, ~, numClusters] = spectralClustering(allPrjs,params{2}); % Pass the number of clusters
                        catch error
                            data = 'Spectral Clustering failed. Aborting. See Console and maybe try again?';
                            fprintf('%s\n',error.message);
                            returnStat = 1;
                            obj.savedState = localStateCopy;
                            return;
                        end
                        % Check we do not go beyond 15 clusters
                        if obj.nClusters + numClusters - numel(selRowsIdx) > 15
                            data = 'Recluster is taking us beyond 15 clusters. Aborting.';
                            returnStat = 1;
                            obj.savedState = localStateCopy;
                            return;
                        end
                        % Assign new IDs
                        if numClusters > numel(selRowsIdx) % More clusters than before
                            newIDs = [obj.displayIDs(selRowsIdx);...
                                ((max(obj.displayIDs)+1):(max(obj.displayIDs) + numClusters - numel(selRowsIdx)))'];
                        else % Less clusters than before
                            newIDs = obj.displayIDs(selRowsIdx(1:numClusters));
                        end
                    else % numel(varargin) == 1, reapplying an old action
                        % Reproduce newIDs, clusterIndexes, numClusters
                        % Then let the magic happen
                        data = varargin{1};
                        newIDs = data{1}(:); % Make sure is a column
                        numClusters = numel(newIDs);
                        clusterIndexes = zeros(numel(allTimes),1);
                        for cl = 1:numClusters
                            [~,r,~] = intersect(allTimes,data{2}{cl});
                            clusterIndexes(r) = cl;
                        end
                    end % varargin switch
                    
                    % Process local variables
                    obj.nClusters = obj.nClusters + numClusters - numel(selRowsIdx);
                    
                    obj.displayIDs(selRowsIdx) = [];
                    obj.displayIDs = [obj.displayIDs ; newIDs];
                    obj.spikeTrains(selRowsIdx) = [];
                    newSpikeTrains = cellfun(@(x) allTimes(clusterIndexes == x),num2cell(1:numClusters)','uni',false);
                    obj.spikeTrains = [obj.spikeTrains; newSpikeTrains];
                    obj.prjTrains(selRowsIdx) = [];
                    newPrjTrains = cellfun(@(x) allPrjs(clusterIndexes == x,:),num2cell(1:numClusters)','uni',false);
                    obj.prjTrains = [obj.prjTrains; newPrjTrains];
                    
                    obj.spikeCounts(selRowsIdx) = [];
                    obj.spikeCounts = [obj.spikeCounts ;...
                        cellfun(@(x) numel(x),newSpikeTrains,'uni',true)];
                    
                    obj.contaminationValues(selRowsIdx) = [];
                    obj.contaminationValues = [obj.contaminationValues;...
                        cellfun(@(x) edu.ucsc.neurobiology.vision.anf.NeuronCleaning.getContam(x,int32(obj.nSamples)),...
                        newSpikeTrains,'uni',true)];
                    
                    obj.statusRaw(selRowsIdx,:) = [];
                    obj.statusRaw = [obj.statusRaw;...
                        repmat([0 0],numClusters,1)];
                    
                    obj.comment(selRowsIdx) = [];
                    obj.comment = [obj.comment ; repmat({'Reclustered'},numClusters,1) ];
                    
                    obj.eisLoaded(selRowsIdx) = [];
                    obj.eisLoaded = [obj.eisLoaded ; cell(numClusters,1)];
                    
                    obj.stasLoaded(selRowsIdx) = [];
                    obj.stasLoaded = [obj.stasLoaded ; cell(numClusters,1)];
                    
                    obj.EIdistMatrix(selRowsIdx,:) = [];
                    obj.EIdistMatrix(:,selRowsIdx) = [];
                    obj.EIdistMatrix = [obj.EIdistMatrix, nan(size(obj.EIdistMatrix,1),numClusters) ;...
                        nan(numClusters,obj.nClusters)];
                    
                    data = {newIDs, newSpikeTrains}; % Add in meaningful data reporting here
                case EditAction.SHRINK % Shrink is cheap, we do not check if we are computing or reapplying
                    % Check ID validity
                    [~,selRowsIdx,~] = intersect(obj.displayIDs,params{1});
                    if numel(selRowsIdx) < numel(params{1})
                        data = 'Some requested IDs are invalid. Aborting edition.';
                        returnStat = 1;
                        obj.savedState = localStateCopy;
                        return;
                    end
                    
                    % Check we do not go beyond 15 clusters
                    if obj.nClusters == 15
                        data = 'Shrinking (inc. outlier cluster) is taking us beyond 15 clusters. Aborting.';
                        returnStat = 1;
                        obj.savedState = localStateCopy;
                        return;
                    end
                    % Process shrink
                    
                    outlierSpikeTrain = zeros(1,0);
                    outlierPrjTrain = zeros(0,5);
                    for c = selRowsIdx'
                        keep = morphShrinkToTarget(obj.spikeTrains{c},obj.prjTrains{c},params{2},params{3},obj.nSamples);
                        outlierSpikeTrain = [outlierSpikeTrain,obj.spikeTrains{c}(~keep)];
                        outlierPrjTrain = [outlierPrjTrain;obj.prjTrains{c}(~keep,:)];
                        obj.spikeTrains{c} = obj.spikeTrains{c}(keep);
                        obj.prjTrains{c} = obj.prjTrains{c}(keep,:);
                    end
                    [outlierSpikeTrain,order] = sort(outlierSpikeTrain);
                    outlierPrjTrain = outlierPrjTrain(order,:);
                    
                    % Process local variables
                    obj.nClusters = obj.nClusters + 1;
                    
                    obj.displayIDs = [obj.displayIDs ; max(obj.displayIDs)+1];
                    obj.spikeTrains = [obj.spikeTrains; {outlierSpikeTrain}];
                    obj.prjTrains = [obj.prjTrains; {outlierPrjTrain}];
                    
                    storeSpikeCounts = obj.spikeCounts(selRowsIdx);
                    obj.spikeCounts(selRowsIdx) = cellfun(@numel, obj.spikeTrains(selRowsIdx));
                    obj.spikeCounts = [obj.spikeCounts ; numel(outlierSpikeTrain)];
                    
                    obj.contaminationValues(selRowsIdx) = cellfun(@(x) edu.ucsc.neurobiology.vision.anf.NeuronCleaning.getContam(x,int32(obj.nSamples)),...
                        obj.spikeTrains(selRowsIdx),'uni',true);
                    obj.contaminationValues = [obj.contaminationValues;...
                        edu.ucsc.neurobiology.vision.anf.NeuronCleaning.getContam(obj.spikeTrains{end},int32(obj.nSamples))];
                    
                    obj.statusRaw(selRowsIdx,:) = 0;
                    obj.statusRaw = [obj.statusRaw;...
                        [-2 , 0]];
                    
                    obj.comment(selRowsIdx) = repmat({'Shrunk down'},numel(selRowsIdx),1);
                    obj.comment = [obj.comment ; {'Shrink outliers - Will be deleted at next consolidate'} ];
                    
                    obj.eisLoaded(selRowsIdx) = {[]};
                    obj.eisLoaded = [obj.eisLoaded ; {[]}];
                    
                    obj.stasLoaded(selRowsIdx) = {[]};
                    obj.stasLoaded = [obj.stasLoaded ; {[]}];
                    
                    obj.EIdistMatrix(selRowsIdx,:) = nan;
                    obj.EIdistMatrix(:,selRowsIdx) = nan;
                    obj.EIdistMatrix = [obj.EIdistMatrix, nan(size(obj.EIdistMatrix,1),1) ;...
                        nan(1,obj.nClusters)];
                    
                    data = {obj.spikeCounts(selRowsIdx)./storeSpikeCounts, obj.contaminationValues(selRowsIdx),...
                        obj.spikeTrains(selRowsIdx)};
                    
                otherwise
                    obj.savedState = localStateCopy;
                    throw(MException('','ClusterEditBackend:localApplyAction - Unhandled EditAction in switch statement.'));
            end
            returnStat = 0;
        end
        
        % Destructor - close handles to java files
        function delete(obj)
            if ~obj.typeIsSpectra
                obj.prj.javafile.close();
                obj.neurons.javafile.close();
            end
            if numel(obj.eiFile) > 0
                obj.eiFile.close();
            end
            if numel(obj.eiFileRaw) > 0
                obj.eiFileRaw.close();
            end
            if numel(obj.staFile) > 0
                obj.staFile.close();
            end
            if numel(obj.staFileRaw) > 0
                obj.staFileRaw.close();
            end
            if numel(obj.globalsFile) > 0
                obj.globalsFile.close();
            end
        end
        
        % function saveLocalState
        %   Handles saving of loaded state into a struct
        %   for purpose of canceling the last action done
        %   Input:
        %       what: String 'clear', 'save' or 'restore'.
        %   Return:
        %       s: 0 if OK, 1 otherwise
        %       msg: Error message for the GUI
        function [s,msg] = saveLocalState(obj,what)
            switch what
                case 'clear'
                    obj.savedState = [];
                case 'save'
                    obj.savedState.nClusters = obj.nClusters;
                    obj.savedState.displayIDs = obj.displayIDs;
                    obj.savedState.contaminationValues = obj.contaminationValues;
                    obj.savedState.spikeCounts = obj.spikeCounts;
                    obj.savedState.statusRaw = obj.statusRaw;
                    obj.savedState.comment = obj.comment;
                    
                    obj.savedState.spikeTrains = obj.spikeTrains;
                    obj.savedState.prjTrains = obj.prjTrains;
                    obj.savedState.eisLoaded = obj.eisLoaded;
                    obj.savedState.EIdistMatrix = obj.EIdistMatrix;
                    obj.savedState.stasLoaded = obj.stasLoaded;
                case 'restore'
                    if numel(obj.savedState) == 0
                        s = 1;
                        msg = 'Nothing to restore.';
                        return;
                    else
                        obj.nClusters = obj.savedState.nClusters;
                        obj.displayIDs = obj.savedState.displayIDs;
                        obj.contaminationValues = obj.savedState.contaminationValues;
                        obj.spikeCounts = obj.savedState.spikeCounts;
                        obj.statusRaw = obj.savedState.statusRaw;
                        obj.comment = obj.savedState.comment;
                        
                        obj.spikeTrains = obj.savedState.spikeTrains;
                        obj.prjTrains = obj.savedState.prjTrains;
                        obj.eisLoaded = obj.savedState.eisLoaded;
                        obj.EIdistMatrix = obj.savedState.EIdistMatrix;
                        obj.stasLoaded = obj.savedState.stasLoaded;
                        
                        obj.savedState = [];
                    end
                otherwise
                    throw(MException('',sprintf('ClusterEditBackend:saveLocalState - Call is not OK with arg %s',what)));
            end
            s = 0;
            msg = '';
        end
    end % End of dynamic methods
    
    
    methods(Static)
        % Statically computes neuron IDs for electrode-clusterID pairs
        % Allows to store absolute neuron information during neuron cleaning
        % Rather than the pattern of successive removals
        % See duplicateRemoval.m, variables saved in cleanPattern.mat
        function IDs = getIDs(el, clust)
            cfg = mVisionConfig();
            cfg = cfg.getSpectralConfig();
            maxClust = cfg.maxEV;
            IDs = (el-2)*maxClust + clust;
        end
        
        function [el, clust] = getElClust(ID)
            cfg = mVisionConfig();
            cfg = cfg.getSpectralConfig();
            maxClust = cfg.maxEV;
            el = floor((ID-1)./maxClust) + 2;
            clust = mod(ID-1,maxClust) + 1;
        end
        
        % function shortenDuplicatePairs
        %   shortens the duplicate discard paths for an array in duplicate removal syntax:
        %   where each row is [ID of neuron 1, ID of neuron 2 discarded as a duplicate of neuron 1]
        %   Shortens the chains meaning that the left columns finally only contains neurons that do
        %   not appear in the second column (ie the only duplicate neuron that is finally kept)
        %
        %   Input:
        %       IDPairs: nx2 array of duplicate discard pairs
        %       rows such that [ID, ID discarded as duplicate]
        function IDPairs = shortenDuplicatesPath(IDPairs)
            % 1st column - ID for which it is deleted
            % 2nd column - ID deleted
            while true
                [~,i,j] = intersect(IDPairs(:,1),IDPairs(:,2));
                if numel(i) == 0
                    break
                end
                IDPairs(i,1) = IDPairs(j,1);
            end
        end
        
    end % Static methods
end % Class
