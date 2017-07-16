% Script to combine signalparams.mat files into an excel spreadsheet
% MS 2016.07.24

names = fieldnames(signalParams)';

if isempty(output)
    output = names;
end

tmpArray = [];
for i = 1:length(names)
   val =  getfield(signalParams,names{i});
   
   tmpArray = [tmpArray,val];
   
end

output = [output;tmpArray];

clearvars -except 'output'