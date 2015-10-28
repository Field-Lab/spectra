function mergePrj( rootFolder, datasets, timeTags )
    %MERGEPRJ merges prj files across datasets
    
    nDatasets = numel(datasets);
    
    info = whos('-file',[datasets{1},'.prj.mat'],'spikeTimes');
    nElectrodes = info.size(1);
    
    nOffset = zeros(1,nDatasets);
    for d = 1:(nDatasets-1)
        load([datasets{d},'.spikes.mat'],'nSamples');
        nOffset(d+1) = nSamples;
    end
    nOffset = cumsum(nOffset);
    nOffset(1) = -str2double(timeTags{1}(2:find(timeTags{1} == '-',1)-1)) * 20000; % HARDCODED SAMPLING RATE
    
    globalSpikes = cell(nElectrodes,1);
    for el = 1:nElectrodes
        globalSpikes{el} = cell(1,nDatasets);
    end
    
    for d = 1:nDatasets
        load([datasets{d},'.prj.mat'],'spikeTimes');
        for el = 1:nElectrodes
            globalSpikes{el}{d} = bsxfun(@plus,spikeTimes{el},nOffset(d));
        end
    end
    spikeTimes = cellfun(@(x) horzcat(x{:}),globalSpikes,'uni', false);
    clear globalSpikes;
    
    load([datasets{1},'.prj.mat'],'eigenValues','eigenVectors');
    save([rootFolder,filesep,'concat.prj.mat'],'eigenValues','eigenVectors','spikeTimes','-v7.3');
    
    % Now working on projections
    for el = 1:nElectrodes
        concatBuffer = cell(1,nDatasets);
        for d = 1:nDataset
            load([datasets{d},'.prj.mat'],sprintf('projSpikes%u',el));
            eval(sprintf('concatBuffer{%u} = projSpikes%u',d,el));
        end
        eval('projSpikes%u = vertcat(concatBuffer{:})',el);
        save([rootFolder,filesep,'concat.prj.mat'],sprintf('projSpikes%u',el),'-append');
    end
end

