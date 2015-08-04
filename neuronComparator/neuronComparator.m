function score = neuronComparator(varargin)
    %NEURONCOMPARATOR loads 2 neurons-raw file and provides the clustering score
    %
    % Loads a reference and a to-compare .neurons-raw files
    % Provides the score of the new clustering
    % versus the reference one by matching sets of clustered spikes
    % Computes matching spikes fraction and maximizes over cluster permutations
    
    javaaddpath([pwd,filesep,'vision',filesep,'Vision.jar']);
    javaaddpath([pwd,filesep,'neuronComparator',filesep]);
    import .*
    
    %%
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
    
    %%
    
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
    
    %%
    maxEl = max(max(neurElRef),max(neurElComp)) + 1;
    score = zeros(1,maxEl);
    
    for el = 2:maxEl
        %%
        seekRef = find(neurElRef == el);
        seekComp = find(neurElComp == el);
        
        if numel(seekRef) == 0 || numel(seekComp) == 0
            continue
        end
        
        setInters = double(neuronCompLoop.computeIntersects(...
            vertcat(neurTimesRef{seekRef}),...
            vertcat(neurTimesComp{seekComp}),...
            neurNumRef(seekRef),...
            neurNumComp(seekComp)));
        
        setInters = bsxfun(@rdivide,setInters,neurNumRef(seekRef));
        % setInters(setInters > 1) = 1./setInters(setInters > 1).^2;
        
        score(el) = neuronCompLoop.maxMetric(setInters)/numel(seekRef);
    end
end
