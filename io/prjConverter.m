%PRJCONVERTER
%
% This script converts .prj.mat files from before 08/10/2015 into the newer format
% Old .prj.mat files: one big nElectrodes x 1 cell array projSpikes, each cell containing spike
% projections
% New .prj.mat files: nElectrodes variables projSpikesXXX containing projections arrays
%
% List of files to convert can be found in shell (find . -name '*.prj.mat', then filtered by date)
% and put in prjMat.input
%
% Note: there should be no old file left
% Note: this script is overwriting - which must be changed -
%   in a manner that if matlab is closed/crashed during execution, all projections of untreated yet
%   electrodes are lost
%
% Author -- Vincent Deo -- Stanford University -- August 21, 2015

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