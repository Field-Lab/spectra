function diaged = neuronComparator(varargin)

%% Script neuron comparator
% Loads 2 neuron file and compares matching fractions of spike times

javaaddpath([pwd,filesep,'vision',filesep,'Vision.jar']);
javaaddpath([pwd,filesep,'neuronComparator',filesep]);
import .*

%%
if nargin == 0
    aName = '';
    neurPath1 = 'X:\EJgroup_data\TestOut\vision_neurons\2005-04-26-0\data002\data002.neurons';
    neurPath2 = 'X:\EJgroup_data\TestOut\matlab_neurons\2005-04-26-0\data002\data002.neurons';
end
if nargin == 1
    aName = [varargin{1}(1:12),filesep,varargin{1}(14:20),filesep,varargin{1}(14:20)];
%     neurPath1 = ['/home/vision/vincent/out_neur_only/',aName,'.neurons'];
%     neurPath2 = ['/home/vision/vincent/ref_neur_only/',aName,'.neurons'];

    neurPath1 = ['X:\EJgroup_data\TestOut\matlab_neurons_3\',aName,'.neurons'];
    neurPath2 = ['X:\EJgroup_data\TestOut\vision_neurons\',aName,'.neurons'];
end
    


neurFile1 = edu.ucsc.neurobiology.vision.io.NeuronFile(neurPath1);
neurFile2 = edu.ucsc.neurobiology.vision.io.NeuronFile(neurPath2);

neurList1 = neurFile1.getIDList();
neurList2 = neurFile2.getIDList();

neurTimes1 = cell(numel(neurList1),1);
neurTimes2 = cell(numel(neurList2),1);

neurNum1 = zeros(numel(neurList1),1);
neurNum2 = zeros(numel(neurList2),1)';


%%
neurCat1 = [];
neurCat2 = [];

for i = 1:numel(neurList1)
    neurTimes1{i} = neurFile1.getSpikeTimes(neurList1(i));
    neurNum1(i) = numel(neurTimes1{i});
    neurCat1 = [neurCat1;neurTimes1{i}];
end
for i = 1:numel(neurList2)
    neurTimes2{i} = neurFile2.getSpikeTimes(neurList2(i));
    neurNum2(i) = numel(neurTimes2{i});
    neurCat2 = [neurCat2;neurTimes2{i}];
end
%%
% setInters = zeros(numel(neurList1),numel(neurList2));
% tic
% for i = 1:numel(neurList1)
%     for j = 1:numel(neurList2)
%         if numel(intersect(neurTimes1{i}(1:99),neurTimes2{j}(1:99))) > 0
%             setInters(i,j) = numel(intersect(neurTimes1{i},neurTimes2{j}));
%         end
%     end
% end
% toc
tic
setInters = double(neuronCompLoop.computeIntersects(neurCat1,neurCat2,neurNum1,neurNum2));
toc

%%
img = cat(3,bsxfun(@rdivide,setInters,neurNum1),repmat(bsxfun(@rdivide,setInters,neurNum2),[1,1,2]));

compare = max(img(:,:,1),img(:,:,2));
%%
[elem,index] = sort(compare(:),'descend');
[row,col] = ind2sub(size(compare),index);

x = unique(row,'stable');
y = unique(col,'stable');

diaged = img(x,y,:);

%%
% figure(numb)
% imshow(diaged);

bName = aName(1:(end-8));
bName(bName == filesep) = '-';
imwrite(diaged,['X:\EJgroup_data\TestOut\comp_neurons_img\',bName,'.png'],'png');
% imwrite(diaged,['/home/vision/vincent/pngCompares/',bName,'.png'],'png');
disp([aName(1:(end-8)),' done.']);
pause(0.1);
end
