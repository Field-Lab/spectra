classdef ClusterEditBackend < handle
    %CLUSTEREDITBACKEND Handles persistently the backend of the cluster visualizer
    
    properties (GetAccess = public, SetAccess = immutable)
        typeIsSpectra
        
        analysisPath
        
        prj % struct - {bool exists ; String path ; matfile or javafile}
        neurons % struct - {bool exists ; String path ; matfile or javafile}
        
        config
        
        % Java related stuff for ei and sta display
        eiFile
        staFile
        globalsFile
        
        arrayID
        electrodeMap
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
        spikeTrainCorr % Matrix of cross correlations between spike trains (actually irrelevant on a single electrode)
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
            else
                obj.typeIsSpectra = 0;
            end
            
            % Type specific loading
            if obj.typeIsSpectra
                % Projections file
                files = dir([analysisPath,filesep,'*.prj.mat']);
                if numel(files) == 1
                    obj.prj.exist = true;
                    obj.prj.path = [analysisPath,filesep,files(1).name];
                else if numel(files) == 0
                        throw(MException('','ClusterEditBackend:ClusterEditBackend - No projections file found'));
                    else
                        throw(MException('','ClusterEditBackend:ClusterEditBackend - Multiple projections files found\nPlease use extended constructor call\nneuronViewer(analysisPath,projectionsFile,neuronsRawFile,modelFile)'));
                    end
                end
                
                % Neurons file
                files = dir([analysisPath,filesep,'*.neurons.mat']);
                if numel(files) == 1
                    obj.neurons.exist = true;
                    obj.neurons.path = [analysisPath,filesep,files(1).name];
                else if numel(files) == 0
                        throw(MException('','ClusterEditBackend:ClusterEditBackend - No neurons-raw file found'));
                    else
                        throw(MException('','ClusterEditBackend:ClusterEditBackend - Multiple neurons-raw files found\nPlease use extended constructor call\nneuronViewer(analysisPath,projectionsFile,neuronsRawFile,modelFile)'));
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
                
                % Clean pattern
                files = dir([analysisPath,filesep,'*.clean.mat']);
                obj.neuronStatuses = zeros(obj.nNeurons,2);
                if numel(files) == 1
                    cleanPatternPath = [analysisPath,filesep,files(1).name];
                    load(cleanPatternPath);
                    IDsDuplicatesRemoved = ClusterEditBackend.shortenDuplicatesPath(IDsDuplicatesRemoved);
                    [~,i,~] = intersect(obj.neuronIDs,IDsRemovedAtContam);
                    obj.neuronStatuses(i,1) = 1; % 1 - removed at contam
                    [~,i,j] = intersect(obj.neuronIDs,IDsMerged(:,2));
                    obj.neuronStatuses(i,1) = 2; % 2 - merged
                    obj.neuronStatuses(i,2) = IDsMerged(j,1);
                    if size(IDsDuplicatesRemoved,2) == 1 % old format, no track
                        [~,i,~] = intersect(obj.neuronIDs,IDsDuplicatesRemoved(:,1));
                        obj.neuronStatuses(i,1) = 3; % 3 - duplicate
                        obj.neuronStatuses(i,2) = -1;
                    else
                        [~,i,j] = intersect(obj.neuronIDs,IDsDuplicatesRemoved(:,2));
                        obj.neuronStatuses(i,1) = 3; % 3 - duplicate
                        obj.neuronStatuses(i,2) = IDsDuplicatesRemoved(j,1);
                    end
                else % no cleaning info found
                    obj.neuronStatuses(:,1) = -1;
                end
            else % Type is vision
                % Projections file
                files = dir([analysisPath,filesep,'*.prj']);
                if numel(files) == 1
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
                    obj.neurons.exist = true;
                    obj.neurons.path = [analysisPath,filesep,files(1).name];
                else if numel(files) == 0
                        throw(MException('','ClusterEditBackend:ClusterEditBackend - No neurons-raw file found'));
                    else
                        throw(MException('','ClusterEditBackend:ClusterEditBackend - Multiple neurons-raw files found\nPlease use extended constructor call\nneuronViewer(analysisPath,projectionsFile,neuronsRawFile,modelFile)'));
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
                    neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile([analysisPath,filesep,files(1).name]);
                    IDsKept = double(neuronFile.getIDList());
                    neuronFile.close();
                    [~,i,~] = intersect(obj.neuronIDs,IDsKept);
                    obj.neuronStatuses(:,1) = -2;
                    obj.neuronStatuses(i,1) = 0;
                else
                    obj.neuronStatuses(:,1) = -1;
                end
            end % type switch
            
            % Initialize EI
            files = dir([analysisPath,filesep,'*.ei']);
            if numel(files) > 0
                eiPath = [analysisPath,filesep,files(1).name];
                obj.eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(eiPath);
            else
                obj.eiFile = []; % Tag for existence of ei File
                fprintf('ClusterEditBackend:ClusterEditBackend - No EI file found, skipping\n');
            end
            
            % Initialize STA
            files = dir([analysisPath,filesep,'*.sta']);
            if numel(files) > 0
                staPath = [analysisPath,filesep,files(1).name];
                obj.staFile = edu.ucsc.neurobiology.vision.io.STAFile(staPath);
            else
                obj.staFile = []; % Tag for existence of ei File
                fprintf('ClusterEditBackend:ClusterEditBackend - No STA file found, skipping\n');
            end
            
            % Grab the globals file
            files = dir([analysisPath,filesep,'*.globals']);
            if numel(files) > 0
                globalsPath = [analysisPath,filesep,files(1).name];
                obj.globalsFile = edu.ucsc.neurobiology.vision.io.chunk.GlobalsFile(globalsPath);
            else
                obj.globalsFile = [];
                obj.staFile = []; % Tag for existence of ei File
                fprintf('ClusterEditBackend:ClusterEditBackend - No Globals file found, STAs Unavailable\n');
            end
            
            % Try to find a classification file
            obj.classification = cell(obj.nNeurons,1);
            files = dir([analysisPath,filesep,'*.classification.txt']);
            if numel(files) == 1
                fid = fopen([analysisPath,filesep,files(1).name]);
                classesRaw = textscan(fid, '%u All/%s', 'delimiter', '\n');
                [~,pos,pos2] = intersect(obj.neuronIDs,classesRaw{1});
                obj.classification(pos) = classesRaw{2}(pos2);
                fclose(fid);
            else
                fprintf('ClusterEditBackend:ClusterEditBackend - Can''t find classification.txt file.\n');
            end
            
            % electrode map
            if numel(obj.eiFile) > 0
                obj.arrayID = obj.eiFile.arrayID;
                obj.electrodeMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(obj.arrayID);
            end
            
            % Handling case where no GUI caller exists - user mode
            % Code will pass through calls to the GUI status bar no questions asked.
            obj.statusBarHandle.String = '';
        end
        
        % function loadEl
        %   loads the requested electrode and prepares all front-end friendly data
        %   Input:
        %       el: MATLAB numbered electrode number
        %
        %   Output:
        %       returnStatus: 1 if load failed (same as already loaded, invalid)
        %                     0 if load succeded
        function returnStatus = loadEl(obj,el)
            if el == obj.elLoaded % Already loaded - save the work
                obj.statusBarHandle.String = sprintf('Electrode %u already loaded.',el-1);
                returnStatus = 1;
                return;
            end
            if (el <= 1) || (el > obj.nElectrodes) % Invalid number
                obj.statusBarHandle.String = sprintf('Electrode number %u invalid, nothing done.',el-1);
                returnStatus = 1;
                return;
            end
            
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
            % Load neuron time correlations
            obj.loadCorrelations();
            
            % Finish
            obj.isDataReady = true;
            returnStatus = 0;
        end
        
        % function loadEI
        %   loads the EIs from the EI file for the IDs in obj.displayIDs
        function loadEI(obj)
            obj.eisLoaded = cell(obj.nClusters,1);
            for c = 1:obj.nClusters
                try
                    obj.eisLoaded{c} = obj.eiFile.getImage(obj.displayIDs(c));
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
        
        % function loadCorrelations
        %   Computes the time correlation between spike trains loaded in obj.spikeTrains
        %   Using parameters corresponding to Vision's duplicate removal defaults
        %   This information is not displayed anywhere at this point, because irrelevant
        %   for neurons on the same electrode
        function loadCorrelations(obj)
            visionCoincidenceTime = 10;
            obj.spikeTrainCorr = ones(obj.nClusters);
            for i = 1:obj.nClusters
                for j = (i+1):obj.nClusters
                    obj.spikeTrainCorr(i,j) = ...
                        edu.ucsc.neurobiology.vision.anf.NeuronCleaning.getCorrVal(...
                        obj.spikeTrains{i},obj.spikeTrains{j},...
                        visionCoincidenceTime);
                    obj.spikeTrainCorr(j,i) = obj.spikeTrainCorr(i,j);
                end
            end
        end
        
        % Function loadSTA
        %   load the STAs in the STA file for the IDs in obj.displayIDs
        function loadSTA(obj)
            obj.stasLoaded = cell(obj.nClusters,1);
            for c = 1:obj.nClusters
                try
                    obj.stasLoaded{c} = obj.staFile.getSTA(obj.displayIDs(c));
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
        
        % function softApplyAction
        %   applies an edit action to the data stored in the cache of this object,
        %   locally, on the selected electrode.
        %   This is previewing feature for edition actions.
        %
        %   Inputs:
        %       action: EditAction object
        %       parameters: cell array describing action parameters
        function softApplyAction(obj,action,parameters)
            validateattributes(action,{'EditAction'},{});
            validateattributes(parameters,{'cell'},{});
        end
        
        % function hardApplyAction
        %   Overwrites the current state of displayed clusters to the neurons.mat file
        %
        %   Inputs:
        %       action: EditAction object
        %       parameters: cell array describing action parameters
        function hardApplyPreviewedActions(obj,action,parameters)
            validateattributes(action,{'EditAction'},{});
            validateattributes(parameters,{'cell'},{});
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
