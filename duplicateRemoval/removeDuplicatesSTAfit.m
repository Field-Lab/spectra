function removeDuplicatesEIComparison(dataFolder, varargin)
    %REMOVEDUPLICATESSTAFIT Summary of this function goes here
    %
    % This function identifies duplicates in a neuron-raw file, based on
    % Bhattacharyya distance between gaussian STA fits.
    % This function then saves a duplicate removed .neurons file in the dataset folder
    % If a neurons file is already present, no overwriting is done unless
    % 'true' is specified as second argument.
    % This requires to compute STAs and parameters file on the neurons-raw file
    % (rename .neurons-raw into .neurons and run appropriate vision calculation,
    % then rename .params into .params-raw-green)
    %
    % Author -- Vincent Deo -- Stanford University -- September 21, 2015
    
    timer = tic;
    hyperThreshold = -log(0.98); % Merge at given overlapping probability.
    hyperMaxAxisRatio = 4;
    
    narginchk(1,3);
    if nargin >= 2
        overwrite = varargin{1};
        validateattributes(overwrite,{'logical'},{},'','overwriteTag',2);
    else
        overwrite = false;
    end
    if nargin == 3
        preclean = varargin{2};
        validateattributes(preclean,{'logical'},{},'','preCleaningTag',2);
    else
        preclean = false;
    end
        
    
    % Checking java path
    javaaddpath ./vision/
    javaaddpath ./duplicateRemoval/java_EI_comparison/
    javaaddpath ./vision/Vision.jar -end
    
    % Making sure dataFolder ends by '\' or '/', whichever is right
    if dataFolder(end:end)~=filesep
        dataFolder = [dataFolder filesep];
    end
    
    % Finding checking files
    list = dir([dataFolder,'*.neurons-raw']);
    if numel(list) == 1
        neurRawPath = list.name;
        datasetName = neurRawPath(1:find(neurRawPath == '.',1,'last')-1);
        neurRawPath = [dataFolder,neurRawPath];
    else
        throw(MException('','None or several neurons-raw found.'));
    end
    
    neurPath = [dataFolder datasetName '.neurons'];
    paramsPath = [dataFolder datasetName '.params-raw-green'];
%     paramsPath = [dataFolder datasetName '.params-cleaned'];
    
    list = dir(neurPath);
    if numel(list) >= 1
        if overwrite
            fprintf('Warning: %s is overwritten.\n',list.name);
            delete(neurPath)
        else
            fprintf('Warning: %s should be overwritten at the end of this execution, but overwrite is disabled. Returning now.\n',list.name);
            return
        end
    end
    
    if preclean % Precleaning is requested - Start low count removal and contamination cut
        system(['java -cp .\vision\Vision.jar edu.ucsc.neurobiology.vision.calculations.CalculationManager -c .\vision\config.xml "Neuron Cleaning" ',neurRawPath,' 100 0.1 0 0.25']);
    end
    
    %% Load parameters file - setting up formats
    % Load Gaussian fits
    paramsFile = edu.ucsc.neurobiology.vision.io.ParametersFile(...
        paramsPath);
    
    % Load neuron IDs
    if preclean
        % Temporary pre-cleaned neuron file
        neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neurPath);
        IDs = neuronFile.getFullIDList();
        neuronFile.close();
        delete(neurPath);
    else
        IDs = paramsFile.getIDList();
    end
    gsn = paramsFile.gatherGaussians(IDs);
%     % gsn comes up in format [x center, y center, x size, y size, theta]
%     % Default config is 5-sig contour. Normalizing to 1-sig
%     gsn(:,[3,4]) = gsn(:,[3,4]) ./ 5;
    discardByAxisRatio = or(gsn(:,4)./gsn(:,3) > hyperMaxAxisRatio,gsn(:,4)./gsn(:,3) < 1./hyperMaxAxisRatio);
    gsn(discardByAxisRatio,:) = [];
    IDs(discardByAxisRatio) = [];
    
    % Plotting utilities, to compare with display in vision
    if false
        %%
        xPlot = @(g,t) bsxfun(@times,cos(g(:,5)).*g(:,3),cos(t)) + bsxfun(@times,sin(g(:,5)).*g(:,4),sin(t));
        yPlot = @(g,t) bsxfun(@times,-sin(g(:,5)).*g(:,3),cos(t)) + bsxfun(@times,cos(g(:,5)).*g(:,4),sin(t));
        samp = 0:0.01:2*pi;
        x = bsxfun(@plus,xPlot(gsn,samp),gsn(:,1));
        y = bsxfun(@plus,yPlot(gsn,samp),gsn(:,2));
        plot(x',y','k','linewidth',0.1);
    end
    % Function to extract ellipses covariance matrices, in n x 4 format, for all rows of gsn.
    covMat = @(g) [g(:,3).^2 .* cos(g(:,5)).^2 + g(:,4).^2 .* sin(g(:,5)).^2,...
        (g(:,3).^2 - g(:,4).^2).*sin(g(:,5)).*cos(g(:,5)) * [1,1],...
        g(:,3).^2 .* sin(g(:,5)).^2 + g(:,4).^2 .* cos(g(:,5)).^2];
    
    % Determinant function, computing for all rows
    det = @(cov) cov(:,1) .* cov(:,4) - cov(:,2) .* cov(:,3);
    
    % Covariance matrices and determinants
    c = covMat(gsn);
    d = det(c);
    
    % Cleaning params file
    paramsFile.close(false); clear paramsFile;
    
    %% Computing Bhattacharyya distance
    % Pairwise average covariance matrice between ellipses
    % The covariance matrix of the pair is stored along the 3rd dimension, in a
    % nNeurons x nNeurons x 4 array;
    pairCov = 0.5 .* bsxfun(@plus,permute(c,[1,3,2]),permute(c,[3,1,2]));
    
    % covariance matrix determinant along 3rd dimension
    det3d = @(cov) cov(:,:,1) .* cov(:,:,4) - cov(:,:,2) .* cov(:,:,3);
    
    % log term of the Bhattacharyya distance
    logTerm = 0.5 .* log(det3d(pairCov)./bsxfun(@times,sqrt(d),sqrt(d)'));
    
    % comatrix inversion of covariance matrices
    pairCov = bsxfun(@rdivide,bsxfun(@times,pairCov(:,:,[1,3,2,4]),permute([1,-1,-1,1],[3,1,2])),det3d(pairCov));
    
    % pairwise difference between ellipses centroids
    muDiffs = bsxfun(@minus,permute(gsn(:,1:2),[1,3,2]),permute(gsn(:,1:2),[3,1,2]));
    
    % 3rd dimension quadratic product computation
    % A is a matrix n x n x 4 of n x n different 2x2 matrices
    % x is a matrix n x n x 2 of n x n different 2x1 array
    % Computes n x n matrix of x' * A * x result for all n x n cases
    thirdDimMul = @(x,A) x(:,:,1).^2 .* A(:,:,1) + (x(:,:,1) .* x(:,:,2)) .* (A(:,:,2) + A(:,:,3)) + x(:,:,2).^2 .* A(:,:,4);
    
    % General term of Bhattacharyya distance
    bhatta = 0.125 .* thirdDimMul(muDiffs,pairCov) + logTerm;
    
    %% Extracting distance values
    bhatta(isnan(bhatta)) = 0;
    [row,col] = find( bhatta ); % remove singularities
    % remove diagonal terms if any
    diagRem = row == col;
    row(diagRem) = [];
    col(diagRem) = [];
    % remove duplicate indices for computing each edge only once
    duplicatePairs = unique(sort([row,col],2),'rows');
    
    % List of distance-valued connected neuron pairs
    % [row, col, edge_value]
    pairList = [double(duplicatePairs),bhatta(sub2ind(size(bhatta),duplicatePairs(:,1),duplicatePairs(:,2)))];
    % Sort by increasing edge value - required for graph/conn component analysis
    pairList = sortrows(pairList,3);
    
    %% Threshold-based connected component analysis
    g = Graph(max(max(pairList(:,1:2))));
    g.debug = true;
    
    % Process all edges
    g.addAllEdges(pairList(:,1)-1,pairList(:,2)-1,pairList(:,3));
    
    % Extract final connected component structure and merge edge values list
    g.removeSingletons();
    t = g.getFinalTree();
    thr = t.thresholdList();
    
    %% Rectangle-style display of graph structure
    if false
        %%
        rects = t.enumRectangles();
        figure(1)
        plot([rects(:,1),rects(:,2),rects(:,2),rects(:,1)]',[rects(:,4),rects(:,4),rects(:,3),rects(:,3)]')
        figure(2)
        plot([0;sort(thr,'ascend')]',(numel(thr)+1):-1:1)
    end
    
    %% Threshold setting and tree splitting
    thrCut = hyperThreshold; % Set threshold
    
    % Get connected component partition at selected threshold
    treeList = t.partitioning(thrCut);
    
    % Neuron-raw file
    neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neurRawPath);
    
    % auxiliary functions - to keep max spiking neuron of each equivalence class
    nSpikes = @(x) neuronFile.getSpikeCount(IDs(x));
    argmax = @(x) find(x == max(x));
    
    % Neurons to keep
    keep = cell2mat(cellfun(@(x) x(argmax(arrayfun(nSpikes,x+1))), treeList, 'uni',false));
    
    %% Erase unkept neurons, write .neurons file in dataset folder
    selectNeuronsFromRawFile(dataFolder,IDs(keep+1));
    
    fprintf('Duplicate removal finished for %s.\nTook %u seconds. %u neurons kept from %u.\n',...
        dataFolder(min(find(dataFolder == '\',3,'last')+1):end),round(toc(timer)),numel(keep),numel(IDs));
    
    if false
        %%
        samp = 0:0.01:2*pi;
        x = bsxfun(@plus,xPlot(gsn(keep+1,:),samp),gsn(keep+1,1));
        y = bsxfun(@plus,yPlot(gsn(keep+1,:),samp),gsn(keep+1,2));
        hold on; plot(x',y','r','linewidth',0.1); hold off
        axis([10 50 5 28]);
    end
    
    system(['java -cp .\vision\Vision.jar edu.ucsc.neurobiology.vision.calculations.CalculationManager -c .\vision\config.xml "Make Parameters File" ',dataFolder,' 4 true true 1 5.0 1 false true 3.0 0 0.02551 -40.0 true 100.0 0.5 false 0 0 false 0 0 false 0 false 0 true']);
end

