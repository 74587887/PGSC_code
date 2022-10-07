clear;

currentFolder = pwd;
addpath(genpath(currentFolder));
addpath(genpath('.'));
file=dir(fullfile('./data','*.mat'));


para =[10,0.01,10^-5];
% the optimal parameter set varies with different datasets



for i = 1:length(file)
    load(file(i).name);
    for i=1:size(data,1)
        dist = max(max(data{i})) - min(min(data{i}));
        m01 = (data{i} - min(min(data{i})))/dist;
        data{i} = 2 * m01 - 1;
    end

    [result,Tim]=FGN(data,labels,para(1),para(2),para(3));   
    result    
    dlmwrite('results.txt',result,'precision', '%.5f','-append','delimiter','\t');

end

