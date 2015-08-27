function score = neuronComparator(varargin)
    %NEURONCOMPARATOR loads 2 .neurons-raw files and performs for each electrode
    % a 1-on-1 neuron matching score
    %
    % Loads a reference and a to-compare .neurons-raw files
    % Provides the score of the new clustering
    % versus the reference one by matching sets of clustered spikes
    %
    % Computes matching spikes fraction between clusters pairwise and
    % maximizes over possible cluster permutations
    %
    % Inputs:
    %   0: Runs in script mode - arguments (file paths) hardcoded below
    %   1: Dataset names: Runs in function mode. Argument is the dataset name in
    %       format 'yyyy-mm-dd-x/data0xx'
    %       Folders in which to seek reference and comparison neurons-raw hardcoded below
    %
    % Author -- Vincent Deo -- Stanford University -- August 27, 2015
    
    % Loaders
    javaaddpath([pwd,filesep,'vision',filesep,'Vision.jar']);
    javaaddpath([pwd,filesep,'neuronComparator',filesep]);
    import .*
    
    % Argument parsing
    if nargin == 0 % Script mode - testing
        %%
        aName = '';
        neurPathRef = 'X:\EJgroup_data\NeurRawTest\vision\data000.neurons-raw';
        neurPathComp = 'X:\EJgroup_data\NeurRawtest\matlab\data000.neurons-raw';
    end
    if nargin == 1 % function mode
        %%
        aName = [varargin{1}(1:12),filesep,varargin{1}(14:20),filesep,varargin{1}(14:20)];
        % aName = 'yyyy-mm-dd-xx/data00x/data00x/' or identical with '\'
        
        % My own file system
        % neurPath1 = ['/home/vision/vincent/out_neur_only/',aName,'.neurons'];
        % neurPath2 = ['/home/vision/vincent/ref_neur_only/',aName,'.neurons'];
        
        % Lab file system
        refFolder = '/Volumes/Lab/Projects/spikesorting/mvision/outputsVision/';
        compFolder = '/Volumes/Lab/Projects/spikesorting/mvision/outputs/';
        neurPathRef = [refFolder,aName,'.neurons-raw'];
        neurPathComp = [compFolder,aName,'.neurons-raw'];
    end
    
    % Get files, ID lists, allocate for spike times
    neurFileRef = edu.ucsc.neurobiology.vision.io.NeuronFile(neurPathRef);
    neurFileComp = edu.ucsc.neurobiology.vision.io.NeuronFile(neurPathComp);
    
    neurListRef = neurFileRef.getIDList();
    neurListComp = neurFileComp.getIDList();
    
    neurTimesRef = cell(numel(neurListRef),1);
    neurTimesComp = cell(numel(neurListComp),1);
    
    neurNumRef = zeros(numel(neurListRef),1);
    neurNumComp = zeros(numel(neurListComp),1)';
    
    neurElRef = zeros(numel(neurListRef),1);
    neurElComp = zeros(numel(neurListComp),1)';
    
    %% Load from files
    for i = 1:numel(neurListRef)
        neurTimesRef{i} = neurFileRef.getSpikeTimes(neurListRef(i));
        neurNumRef(i) = numel(neurTimesRef{i});
        neurElRef(i) = neurFileRef.getElectrode(neurListRef(i));
    end
    for i = 1:numel(neurListComp)
        neurTimesComp{i} = neurFileComp.getSpikeTimes(neurListComp(i));
        neurNumComp(i) = numel(neurTimesComp{i});
        neurElComp(i) = neurFileComp.getElectrode(neurListComp(i));
    end
    
    %% Allocate output
    maxEl = max(max(neurElRef),max(neurElComp)) + 1;
    score = zeros(1,maxEl);
    
    %% Comparison loop
    for el = 1:maxEl
        %% Find electrode's neurons (= neuron extraction in vision's terms)
        seekRef = find(neurElRef == el);
        seekComp = find(neurElComp == el);
        
        if numel(seekRef) == 0 || numel(seekComp) == 0
            continue
        end
        
        % Compute pairwise intersects (java)
        setInters = double(neuronCompLoop.computeIntersects(...
            vertcat(neurTimesRef{seekRef}),...
            vertcat(neurTimesComp{seekComp}),...
            neurNumRef(seekRef),...
            neurNumComp(seekComp)));
        
        % Normalize intersects:
        % vs. reference cluster size
        setInters = bsxfun(@rdivide,setInters,neurNumRef(seekRef));
        % vs. reference and new cluster size
        % setInters = sqrt(bsxfun(@rdivide,setInters,neurNumRef(seekRef)).*bsxfun(@rdivide,setInters,neurNumComp(seekComp)));
        
        % Compute best-matching permutation of clusters (java)
        score(el) = neuronCompLoop.maxMetric(setInters)/numel(seekRef);
    end
end
