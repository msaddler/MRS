% startup_MS.m
% MS 2016.06.14

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Setup paths for handling dtag data   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('G:\Mark\code') % <-- Set current directory to MS's code

%Work with blue whale data (bm)
dataset = 'bm15';

% Find right drive
if exist(['E:\' dataset])
    drive='E';
elseif exist(['F:\' dataset])
    drive='F';
elseif exist(['G:\' dataset])
    drive='G';
elseif exist(['H:\' dataset])
    drive='H';
elseif exist(['I:\' dataset])
    drive='I';    
elseif exist(['J:\' dataset])
    drive='J';  
elseif exist(['K:\' dataset])
    drive='K';
else
    drive=[];
    disp('No tagtools installed - check drive or change dataset')
end

if ~isempty(drive)
	% Add path so MatLab can find Tag Tools
	addpath ([drive ':\Mark\dtagtools\d2matlab'])
	addpath ([drive ':\Mark\dtagtools\d3analysis'])
	addpath ([drive ':\Mark\dtagtools\d3matlab_2014'])
    addpath ([drive ':\Mark\dtagtools\x3toolbox'])
    addpath ([drive ':\Mark\dtagtools\x3toolbox\XML4MATv2'])
    addpath ([drive ':\Mark\dtagtools\d3cal'])
    
	% Set paths for DTAG files
	settagpath('cal',[drive ':\Mark\data_bm\cal']);
	settagpath('raw',[drive ':\Mark\data_bm\raw']);
	settagpath('prh',[drive ':\Mark\data_bm\prh']);
	settagpath('audit',[drive ':\Mark\data_bm\audit']) ;
	settagpath('audio',[drive ':\']);
	
	disp(['Tag tools for dataset ' dataset ' running from drive ' drive])
end

%% Quick run dtagaudit

tag = 'bm15_054a';
R = loadaudit_MS(tag);
R = dtagaudit_MS(tag,1,R);