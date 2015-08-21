
fileList = textread(['..',filesep,'prjMatList.input'],'%s');
n = numel(fileList);

for f = 1:n
    fileList{f}
    load(fileList{f});
    if exist('projSpikes','var')
        save(fileList{f},'eigenValues','eigenVectors','spikeTimes','-v7.3');
        for el = 1:numel(projSpikes)
            eval(sprintf('projSpikes%u = projSpikes{el};',el));
            save(fileList{f},sprintf('projSpikes%u',el),'-append');
            eval(sprintf('projSpikes%u = [];',el));
            projSpikes{el} = []; % Progressive RAM clean-up
        end
    end
end