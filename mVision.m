function mVision(dataPath, saveFolder, timeCommand, movieXML, tryToDo, force, varargin)
    %MVISION Main function for single dataset white noise analysis
    %
    % This script manages the successive steps of serial neuron finding
    % Input:
    %   dataSetName: path to the dataset folder
    %       ex: something ending in ...'yyyy-mm-dd-n/dataxxx'
    
    %   saveFolder: the output folder of the analysis.
    %       It is MOST recommended (ie do it) the deepest level subfolder has the same
    %       name 'dataxxx' as the dataset folder. Vision infers file names from
    %       destination folder whereas mVision does it from the source folder.
    %
    %   timeCommand: vision style time command for selection within a dataset
    %       Use empty string '' for full dataset
    %       ex: '', '(10-)', '(500-1000)','(-100)'
    %       Warning: mVision.m DOES NOT support concatenated dataset syntaxes
    %       such as data000(1000-)-data002(-100).
    %       Use concatenatedAnalysis.m (WIP) for that purpose
    %   
    %   movieXML: path to movie XML file
    %       
    %   tryToDo: a 1xnSteps binary array for what subcalculations to do (see list below)
    %       0 - do not try the calculation
    %       1 - try the calculation
    %       String argument 'all' is allowed as shorthand for [1 1 1 1 1 1]
    %       If a calculation is required but the previous was not and the input
    %       files are missing, mVision will crash.
    %
    %   force: a 1xnSteps binary array for overwrite authorization of the calculations
    %       0 - do not overwrite, skip the calculation if its output files are found
    %       1 - overwrite, do the calculation in all cases
    %       String arguments 'all' and 'none' are allowed as shorthand for
    %       respectively all ones and all zeros.
    %
    %   OPTIONALLY
    %   varargin{1} : configTag - a valid configuration tag (string) in the STATIC_CONFIG_LIST file
    %       Otherwise, configuration initializes to defaults as listed in mVisionConfig.m
    %
    %
    % This function successively calls if requested or necessary the calculations:
    %       1/ Noise evaluation
    %       2/ Spike Finding
    %       3/ Covariance Calculation
    %       4/ Projections Calculations
    %       5/ Clustering
    %       6/ Computation of raw EIs and STAs
    %       7/ Automated duplicate removal with stack saving
    %       8/ Export of a .neurons after applying automatic + user cleaning/edit stack
    %       
    % All intermediate files are stored in convenient .mat files.
    % The final .neurons.mat file is converted to a vision compatible .neurons
    %
    % VISION POST PROCESSING
    % To run further analysis with vision, first call the calculations
    % "Generate globals file" then "Copy Raw Data Header To Globals"
    % You can the compute EIs, STAs and the parameter file.
    %
    % Author -- Vincent Deo -- Stanford University -- February 19, 2016
    
    %% VERSION NUMBER %%
    % Increment at each master merge
    version = '0.2';
    fprintf('---------------- Welcome to Spectra v. %s ----------------\n\n',version);
 
    %% SETUP
    % Add subfolders to the matlab path
    addpath(genpath(['.',filesep]));
    
    
    % Create global configuration handle
    narginchk(6,7);
    global GLOBAL_CONFIG
    if numel(varargin) == 1
        GLOBAL_CONFIG = mVisionConfig(varargin{1});
    else % == 0
        GLOBAL_CONFIG = mVisionConfig();
    end
    
    % generate repository path
    repoPath = pwd;
    
    % Java setup - using vision jar available in the repository
    visionPath = [repoPath,filesep,'vision',filesep,'Vision.jar'];
    if ~exist('edu/ucsc/neurobiology/vision/io/RawDataFile','class')
        javaaddpath(visionPath,'end')
    end
    javaaddpath(['.',filesep,'vision']);
    javaaddpath(['.',filesep,'clusterUtil']);
    javaaddpath(['.',filesep,'duplicateRemoval',filesep,'java_EI_comparison',filesep]);
    
    % USER INPUT - Set up data and output folders
    if dataPath(end) == filesep
        dataPath = dataPath(1:(end-1));
    end
    if saveFolder(end) == filesep
        saveFolder = saveFolder(1:(end-1));
    end
    if ~strcmp(dataPath,[filesep,'concat']) && ...
        ~(exist(dataPath,'file') == 2 || exist(dataPath,'file') == 7)
        % [filesep, 'concat'] is a special token here for concatenated analysis
        throw(MException('','demoScript: data source folder|file does not exist'));
    end
    
    mkdir(saveFolder);
    [~,datasetName,~] = fileparts(dataPath); % Catching dataset name as last part of saveFolder
    
    if exist(saveFolder,'file') == 2
        [~,datasetName,~] = fileparts(saveFolder);
    end
    
    % USER input - tryToDo -- won't do any task unless stated here
    nSteps = 8;
    
    % Steps in order:
    % --------- noise - spike - cov - prj - clust - save ----------------------
    if isa(tryToDo,'char')
        tryToDo = lower(tryToDo);
        if strcmp(tryToDo,'all')
            tryToDo = ones(1,nSteps);
        end
        if isa(tryToDo,'char')
            throw(MException('','Invalid tryToDo argument'));
        end
    else
        validateattributes(tryToDo,{'numeric'},{'row','ncols',nSteps,'binary'},'','tryToDo',4);
    end
    
    if isa(force,'char')
        force = lower(force);
        if strcmp(force,'all')
            force = ones(1,nSteps);
        end
        if strcmp(force,'none')
            force = zeros(1,nSteps);
        end
        if isa(force,'char')
            throw(MException('','Invalid force argument'));
        end
    else
        validateattributes(force,{'numeric'},{'row','ncols',nSteps,'binary'},'','force',4);
    end
    
    
    % DEBUG - additional saved file dataset name extension
    nameExt = '';
    
    totalTime = tic;
    
    %% %%%%%%%% Processing steps %%%%%%%%%%
    %
    % Step Header (in case of pipeline modification) - which should also include documentation modif
    % as tryToDo and force arguments change size - should be of the form:
    % thisStep = __STEP_NUMBER__;
    % if tryToDo(thisStep) &&...
    %         (force(thisStep) || ~(exist([saveFolder,filesep,datasetName,'__REQUIRED_OUTPUT_EXTENSION_'],'file') == 2))
    % Also don't forget to increment nSteps somewhere above here.
    
    %% Process noise and make a .noise file
    thisStep = 1;
    if tryToDo(thisStep) &&...
            (force(thisStep) || ~(exist([saveFolder,filesep,datasetName,'.noise'],'file') == 2))
        %%
        fprintf('Starting noise finding...\n');
        tic
        noise = RawDataNoiseEvaluationM(dataPath, saveFolder);
        
        fprintf('Time for noise evaluation %.2f seconds.\n', toc);
    else
        fprintf('Noise not requested or .noise file found - skipping raw data noise evaluation.\n');
    end
    
    
    %% Find spikes and make a .spikes file
    thisStep = 2;
    if tryToDo(thisStep) &&...
            (force(thisStep) || ~(exist([saveFolder,filesep,datasetName,'.spikes.mat'],'file') == 2))
        %%
        fprintf('Starting spike finding...\n');
        tic
        
        sigmaFileName = [saveFolder,filesep,datasetName,'.noise'];
        
        [spikes,ttlTimes,nSamples] = SpikeFindingM(dataPath, saveFolder, timeCommand, sigmaFileName);
        spikeSave = int32(spikes(:,1:2));
        save([saveFolder,filesep,datasetName,'.spikes.mat'],'spikeSave','ttlTimes','nSamples','-v7.3');
        
        fprintf('Time for spike finding %.2f seconds\n',toc);
    else
        fprintf('Spike not requested or .spikes.mat file found - skipping spike finding.\n');
    end
    
    
    %% Covariance calculation
    thisStep = 3;
    if tryToDo(thisStep) &&...
            (force(thisStep) || ~(exist([saveFolder,filesep,datasetName,'.cov.mat'],'file') == 2))
        %%
        fprintf('Starting covariance calculation...\n');
        tic
        
        if ~exist('spikeSave')
            load([saveFolder,filesep,datasetName,'.spikes.mat']);
        end
        
        [covMatrix,averages,totSpikes] = buildCovariances(double(spikeSave), dataPath, timeCommand);
        
        covConfig = GLOBAL_CONFIG.getCovConfig();
        if covConfig.whitening
            noiseSpikes = generateNoiseEvents(double(spikeSave));
           [noiseCovMatrix,~,~] = buildCovariances(noiseSpikes,dataPath,timeCommand);
           covMatrix = whitenMatrices(dataPath, totSpikes, covMatrix, noiseCovMatrix);
%          covMatrix = noiseCovMatrix;
        end
        save([saveFolder,filesep,datasetName,'.cov.mat'],'covMatrix','averages','totSpikes','-v7.3');
        
        fprintf('Time for covariance calculation %.2f seconds\n',toc);
    else
        fprintf('Cov not requested or .cov.mat file found - skipping covariance calculation.\n');
    end
    
    %% Eigenspikes Projections calculation
    thisStep = 4;
    if tryToDo(thisStep) &&...
            (force(thisStep) || ~(exist([saveFolder,filesep,datasetName,'.prj.mat'],'file') == 2))
        %%
        fprintf('Starting projections calculation...\n');
        tic
        if ~exist('covMatrix')
            load([saveFolder,filesep,datasetName,'.cov.mat']);
        end
        
        if ~exist('spikeSave')
            load([saveFolder,filesep,datasetName,'.spikes.mat']);
        end
        
        [projSpikes,eigenValues,eigenVectors,spikeTimes] = ...
            PCProj(dataPath, timeCommand, ...
            double(spikeSave), covMatrix, averages, totSpikes);
        
        save([saveFolder,filesep,datasetName,'.prj.mat'],'eigenValues','eigenVectors','spikeTimes','-v7.3');
        for el = 1:numel(projSpikes)
            eval(sprintf('projSpikes%u = projSpikes{el};',el));
            save([saveFolder,filesep,datasetName,'.prj.mat'],sprintf('projSpikes%u',el),'-append');
            eval(sprintf('projSpikes%u = [];',el));
            projSpikes{el} = []; % Progressive RAM clean-up
        end
        
        fprintf('Time for projections calculation %.2f seconds\n',toc);
    else
        fprintf('Prj not requested or .prj.mat file found - skipping projections calculation.\n');
    end
    
    
    %% Clustering
    thisStep = 5;
    if tryToDo(thisStep) &&...
            (force(thisStep) || ~(exist([saveFolder,filesep,datasetName,'.model.mat'],'file') == 2 &&...
            exist([saveFolder,filesep,datasetName,'.neurons.mat'],'file') == 2))
        %%
        fprintf('Starting clustering...\n')
        tic
        % No projections loader at this level
        % Loader of spikeTimes however
        if ~exist('spikeTimes','var')
            load([saveFolder,filesep,datasetName,'.prj.mat'],'spikeTimes');
        end
        % Clustering framework
        [clusterParams,neuronEls,neuronClusters,neuronSpikeTimes] =...
            PCClustering([saveFolder,filesep,datasetName,'.prj.mat'],spikeTimes);
        
        % Save a model file
        save([saveFolder,filesep,datasetName,'.model.mat'],'clusterParams','-v7.3');
        
        % Save a neurons-raw file (.neurons.mat)
        if ~exist('nSamples')
            load([saveFolder,filesep,datasetName,'.spikes.mat'],'nSamples');
        end
        save([saveFolder,filesep,datasetName,'.neurons.mat'],'neuronEls','neuronClusters','neuronSpikeTimes','nSamples','-v7.3');
        
        fprintf('Time for clustering %.2f seconds\n',toc);
    else
        fprintf('Clust not requested or .neurons|model.mat files found - skipping clustering.\n');
    end
    
    %% Save a neurons raw, compute a *.ei-raw and a *.sta-raw
    thisStep = 6;
    if tryToDo(thisStep)
        computeCfg = GLOBAL_CONFIG.getComputeConfig();
        hasEI = exist([saveFolder,filesep,datasetName,'.ei-raw'],'file') == 2;
        hasSTA = exist([saveFolder,filesep,datasetName,'.sta-raw'],'file') == 2;
        if (computeCfg.ei && ~hasEI) || (computeCfg.sta && ~hasSTA) || force(thisStep)
            %%
            fprintf('Computing raw EIs and STAs...\n');
            tf = {'false','true'};
            fprintf('Configuration is set to: | EIs - %s | STAs - %s\n',tf{computeCfg.ei + 1},tf{computeCfg.sta + 1});
            fprintf('Will overwrite existing *.ei-raw and *.sta-raw if any.\n');
            tic
            % Reload the neurons.mat information if not there
            if ~exist('neuronSpikeTimes','var')
                load([saveFolder,filesep,datasetName,'.neurons.mat']);
            end
            % Build a neurons file
            neuronSaver = NeuronSaverM(dataPath,saveFolder,datasetName,'',0);
            neuronSaver.pushAllNeurons(neuronEls, neuronClusters, neuronSpikeTimes);
            neuronSaver.close();
            
            % Computations. We use system calls for faster computing
            % EIs
            if computeCfg.eiRaw
                commands = GLOBAL_CONFIG.stringifyEICommand(dataPath,saveFolder,datasetName);
                system(commands{1});
                movefile([saveFolder,filesep,datasetName,'.ei'],...
                    [saveFolder,filesep,datasetName,'.ei-raw'])
            end
            % STAs
            if computeCfg.staRaw
                commands = GLOBAL_CONFIG.stringifySTACommand(dataPath, saveFolder, datasetName);
                % Compute "Make White Noise Movie" and "Calculate Auxiliary Parameters"
                % In matlab-JVM mode.
                system(commands{1}); % Generate globals file
                system(commands{2}); % Copy raw data header to globals
                xmlConfig = edu.ucsc.neurobiology.vision.Config(movieXML);
                edu.ucsc.neurobiology.vision.tasks.RunScript.createWhiteNoiseMovie(xmlConfig, saveFolder);
                edu.ucsc.neurobiology.vision.tasks.RunScript.calcAuxParams(xmlConfig, saveFolder);
                system(commands{3}); % STA Calculation Parallel
                pause(5);
                movefile([saveFolder,filesep,datasetName,'.sta'],...
                    [saveFolder,filesep,datasetName,'.sta-raw'])
            end
            % Remove the neurons file
            pause(5);
            if computeCfg.nRaw
                movefile([saveFolder,filesep,datasetName,'.neurons'],...
                    [saveFolder,filesep,datasetName,'.neurons-raw']);
            else
                delete([saveFolder,filesep,datasetName,'.neurons']);
            end
            fprintf('Raw EIs and STAs computation done.\n');
            fprintf('Time for computation %.2f seconds\n', toc);
        else
            fprintf('Raw EIs|STAs not requested or configuration required files found - skipping raw EIs|STAs.\n');
        end
    end
    
    %% Cleaning and saving neurons
    thisStep = 7;
    if tryToDo(thisStep) && ...
            (force(thisStep) || ~(exist([saveFolder,filesep,datasetName,'.clean.mat'],'file') == 2))
        %%
        fprintf('Computing automatic duplicate removal merges/discards...\n')
        tic
        if ~exist('neuronSpikeTimes','var')
            load([saveFolder,filesep,datasetName,'.neurons.mat']);
        end
        
        duplicateRemoval(dataPath, saveFolder, datasetName, ...
            neuronEls, neuronClusters, neuronSpikeTimes);
        
        fprintf('Neuron cleaning done.\n');
        
        fprintf('Time for automatic duplicate removal calculation %.2f seconds\n', toc);
    else
        fprintf('Cleaning pattern not requested or .clean.mat file found - skipping cleaning|saving.\n');
    end
    
    %% Consolidate to latest .neurons.mat/.neurons file, compute eis and stas.
    thisStep = 8;
    if tryToDo(thisStep) && ...
            (force(thisStep) || ~(exist([saveFolder,filesep,datasetName,'.neurons'],'file') == 2)) 
        %%
        fprintf('Exporting automatic and manual edits into a final .neurons file...\nComputing final EIs and STAs...\n')
        tic
        % Load latest neurons
        if exist([saveFolder,filesep,datasetName,'-edited.neurons.mat'],'file') == 2 % Not the first set of manual actions
            fprintf('    Starting from consolidated neurons file.\n');
            load([saveFolder,filesep,datasetName,'-edited.neurons.mat']);
        else
            fprintf('    Starting from raw neurons file.\n');
            if ~exist('neuronSpikeTimes','var')
                load([saveFolder,filesep,datasetName,'.neurons.mat']);
            end
            elevatedStatus = false(size(neuronEls,1),1);
        end
        
        % Check for existing actions or new actions
        if exist([saveFolder,filesep,datasetName,'.edit.mat'],'file') == 2
            fprintf('    A manual edition file .edit.mat was found\n');
            load([saveFolder,filesep,datasetName,'.edit.mat']); % Brings in 'manualActions'
            lastConsolidate = find(cellfun(@(x) x  == EditAction.CONSOLIDATE,manualActions(:,1),'uni',true),1,'last');
            if numel(lastConsolidate) == 0
                lastConsolidate = 0;
            end
            manualActionsNew = manualActions((lastConsolidate+1):end,:);
        else
            fprintf('    No manual edition file .edit.mat was found\n');
            manualActions = cell(0,3);
            manualActionsNew = cell(0,3);
        end
        
        % Case: consolidate all actions
        if size(manualActions,1) > 0 && ...
                ~(exist([saveFolder,filesep,datasetName,'-edited.neurons.mat'],'file') == 2)
            
            fprintf('    No incremental raw neurons file. Building one with all actions.\n');
            [neuronEls, neuronClusters, neuronSpikeTimes, elevatedStatus] = ...
                applyManualCleaningPattern(neuronEls, neuronClusters,...
                neuronSpikeTimes, manualActions, elevatedStatus);
            save([saveFolder,filesep,datasetName,'-edited.neurons.mat'],...
                'neuronEls','neuronClusters','neuronSpikeTimes','nSamples','elevatedStatus','-v7.3');
            
            manualActions = [manualActions ;...
                {EditAction.CONSOLIDATE, {}, {} }];
            save([saveFolder,filesep,datasetName,'.edit.mat'],'manualActions','-v7.3');
        
        % Case: consolidation of only new actions 
        elseif size(manualActionsNew,1) > 0 && ...
                exist([saveFolder,filesep,datasetName,'-edited.neurons.mat'],'file') == 2
            fprintf('    An incremental raw neurons file was found. Applying new actions.\n');
            [neuronEls, neuronClusters, neuronSpikeTimes, elevatedStatus] = ...
                applyManualCleaningPattern(neuronEls, neuronClusters,...
                neuronSpikeTimes, manualActionsNew, elevatedStatus);
            save([saveFolder,filesep,datasetName,'-edited.neurons.mat'],...
                'neuronEls','neuronClusters','neuronSpikeTimes','nSamples','elevatedStatus','-v7.3');
            
            manualActions = [manualActions ;...
                {EditAction.CONSOLIDATE, {}, {} }];
            save([saveFolder,filesep,datasetName,'.edit.mat'],'manualActions','-v7.3');
            
        else % Case: skip to automatic cleaning
            fprintf('    No (new) actions to apply, skipping to applying automatic pattern.\n');
        end
        
        % Automatic actions
        if ~exist('autoActions','var')
            load([saveFolder,filesep,datasetName,'.clean.mat']);
        end
        [neuronEls, neuronClusters, neuronSpikeTimes] = ...
            applyAutoCleaningPattern(neuronEls, neuronClusters, neuronSpikeTimes, autoActions, elevatedStatus);
        
        neuronSaver = NeuronSaverM(dataPath,saveFolder,datasetName,'',0);
        neuronSaver.pushAllNeurons(neuronEls, neuronClusters, neuronSpikeTimes);
        neuronSaver.close();
        
        % Computations. We use system calls for faster computing
        computeCfg = GLOBAL_CONFIG.getComputeConfig();
        % EIs
        if computeCfg.ei
            commands = GLOBAL_CONFIG.stringifyEICommand(dataPath,saveFolder,datasetName);
            system(commands{1});
        end
        % STAs
        if computeCfg.sta
            commands = GLOBAL_CONFIG.stringifySTACommand(dataPath, saveFolder, datasetName);
            % Compute "Make White Noise Movie" and "Calculate Auxiliary Parameters"
            % In matlab-JVM mode.
            if ~exist([saveFolder,filesep,datasetName,'.globals'],'file') == 2
                system(commands{1}); % Generate globals file
                system(commands{2}); % Copy raw data header to globals
            end
            % Evalc for catching text outputs.
            xmlConfig = edu.ucsc.neurobiology.vision.Config(movieXML);
            edu.ucsc.neurobiology.vision.tasks.RunScript.createWhiteNoiseMovie(xmlConfig, saveFolder);
            edu.ucsc.neurobiology.vision.tasks.RunScript.calcAuxParams(xmlConfig, saveFolder);
            delete([saveFolder,filesep,datasetName,'.sta']); % Required
            system(commands{3}); % STA Calculation Parallel
        end
        
        if computeCfg.params
            commands = GLOBAL_CONFIG.stringifyParamsCommand(dataPath, saveFolder, datasetName);
            system(commands{1});
            delete([saveFolder,filesep,datasetName,'.globals']);
            system(commands{2});
            system(commands{3});
        end
            
        fprintf('Incremental edition, neurons, EIs, STAs done.\n');
        
        fprintf('Time for exporting %.2f seconds\n', toc);
    else
        fprintf('Export not requested or .neurons file found - skipping export.\n');
    end
    
    %%
    fprintf('\n');
    fprintf('Total pipeline time %.2f seconds\n',toc(totalTime));
    fprintf([strrep(dataPath,'\','\\'),timeCommand,' finished\n']);
    fprintf('----------------------------------------------------------\n');
    
end
