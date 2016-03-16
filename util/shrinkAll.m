function [ neuronSpikeTimes ] = shrinkAll( neuronEls, neuronSpikeTimes, prjPath )
    %SHRINKALL Summary of this function goes here
    load(prjPath,'spikeTimes');
    for el = 1:size(spikeTimes,1)
        cPos = find(neuronEls == el);
        if numel(cPos) == 0
            continue;
        end
        eval(sprintf('load(%s,''projSpikes%u'');',prjPath,el));
        for pos = cPos(:)'
            [~,spIdx,~] = intersect(spikeTimes{el},neuronSpikeTimes{pos});
            eval(sprintf('prj = projSpikes%u(spIdx,:);',el));
            neuronSpikeTimes{pos} = ...
                neuronSpikeTimes{pos}(morphShrinkToTarget(neuronSpikeTimes{pos},prj,1,0.1,36e6));
        end
        spikeTimes{el} = [];
    end
end

