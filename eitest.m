% EI testing stuff
%%
folder = 'X:\EJgroup_data\whiteCompTest\with\data002';
name = 'data002';

% ei Access
eiPath = sprintf('%s%s%s.ei',folder,filesep,name);
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(eiPath);

allIDs = eiFile.getIDList();
nNeurons = numel(allIDs);
% Load Configuration
cfg = mVisionConfig();
cleanConfig = cfg.getCleanConfig();

% Misc for dup removal algorithm
eiSize = cleanConfig.EILP + cleanConfig.EIRP + 1;
minWindowSize = cleanConfig.globMinWin(2) - cleanConfig.globMinWin(1);
samplingVal = cleanConfig.globMinWin(1):cleanConfig.resamplePitch:cleanConfig.globMinWin(2);
basicValues = 1:(eiSize - minWindowSize);

nElectrodes = eiFile.nElectrodes;
[adjacent1,~] = catchAdjWJava( eiFile, 1 );
[adjacent2,~] = catchAdjWJava( eiFile, 2 );

% EI loader
eiStorage = cell(nNeurons,1); % full remaining neurons EIs
varStorage = cell(nNeurons,1);
isRealigned = false(nNeurons,nElectrodes); % tag to upsample any ID/electrode only once

% Preload all EIs
for n = 1:nNeurons
    eiTemp = eiFile.getImage(allIDs(n));
    eiStorage{n} = squeeze(eiTemp(1,:,:))';
    varStorage{n} = squeeze(eiTemp(2,:,:))';
end
fprintf('EIs loaded.\n');

%% Realignment procedure
for n = 1:nNeurons
    fprintf('%u\n',n);
    for el = 2:nElectrodes;
        interpolant = griddedInterpolant(1:eiSize,eiStorage{n}(:,el),'spline');
        [~, offset] = min(interpolant(samplingVal));
        eiStorage{n}(basicValues,el) =...
            interpolant(basicValues + (offset - 1) * cleanConfig.resamplePitch);
        
        interpolantVar = griddedInterpolant(1:eiSize,varStorage{n}(:,el),'spline');
        varStorage{n}(basicValues,el) = ...
            interpolantVar(basicValues + (offset - 1) * cleanConfig.resamplePitch);
    end
    eiStorage{n} = eiStorage{n}(basicValues,:);
    varStorage{n} = varStorage{n}(basicValues,:);
    
    %eiStorage{n} = eiStorage{n}(:);
    %norm2 = sum(eiStorage{n}.^2);
    
    %eiStorage{n} = eiStorage{n}./sqrt(norm2);
    %varStorage{n} = varStorage{n}(:);
    %varStorage{n} = varStorage{n}./norm2;
    
end
%%
%eiStorage2 = horzcat(eiStorage{:});
%varStorage2 = horzcat(varStorage{:});

[adjacent1,~] = catchAdjWJava( eiFile, 2 );
adjacent1mat = cellfun(@(x) [repmat(x(1),19-numel(x),1);x]'-1,adjacent1,'uni',false);
adjacent1mat = cell2mat(adjacent1mat);
neuronEls = floor(allIDs - 1)/15 + 2;
adjacentsNeuron = adjacent1mat(neuronEls,:);

eiStorage = cat(3,eiStorage{:});
eiStorage = permute(eiStorage,[3 2 1]);
varStorage = cat(3,varStorage{:});
varStorage = permute(varStorage,[3 2 1]);