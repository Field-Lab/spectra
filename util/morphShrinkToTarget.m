function keep = morphShrinkToTarget(time,prj,targetFrac,targetContam,nSamples)
    targetNum = floor(targetFrac * numel(time));
    
    [~,D] = knnsearch(prj,prj,'K',5); % Defaults to 5;
    D = D(:,end);
    [D,i] = sort(D);
    
    % Reduce to only the desired number of spikes
    D = D(1:targetNum);
    i = i(1:targetNum);
    
    % Check contamination value - find maximal size by dichotomy
    
    contam = @(n) edu.ucsc.neurobiology.vision.anf.NeuronCleaning.getContam(sort(time(i(1:n))),int32(nSamples));
    if contam(targetNum) > targetContam
        min = 2;
        max = targetNum;
        val = round(0.5*(min+max));
        while (max-min) > 2
            tmp = val;
            if contam(val) < targetContam
                val = floor(0.5*(val+max));
                min = tmp;
            else
                val = floor(0.5*(val + min));
                max = tmp;
            end
        end
    else
        val = targetNum;
    end
    keep = false(val,1);
    keep(i(1:val)) = true;
end