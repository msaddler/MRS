% Script for plotting SEP2 call spectrograms and waveforms
% MS 2016.08.02

% Load-in savedSignals for desired audio signal

cd('G:\Mark\code')
addpath('G:\Mark\code\tools')

pt = 1024;
ovl = 0.99;

FLIM = [5 500];
YLIM = [-.04 .04];
YTICKS = [-.02 0 .02];

% Custom (non-linear) map
MAP = [1:-.01:0.50, 0.49:-0.03:0 ]';
MAP = [MAP,MAP,MAP];

fs = savedSignals{1,1}.fs;

y = savedSignals{1,1}.y;
y = y - mean(y);
t = savedSignals{1,1}.T;

%close all
f = figure('Position',[50 100 1500 600]);
wave_ax = axes('Parent',f,'Position',[.10 .85 .85 .10]);
spec_ax = axes('Parent',f,'Position',[.10 .55 .85 .30]);

MS_spectrogram(spec_ax,y,fs,pt,ovl,FLIM,MAP);
MS_waveform(wave_ax,t,y,YLIM);

linkaxes([spec_ax, wave_ax],'x')

% Clear uneccessary x-axis labels
set(wave_ax,'xticklabel',[])
set(wave_ax,'ytick',YTICKS)

ylabel(spec_ax,{'Frequency';'(Hz)'},'FontSize',10)
ylabel(wave_ax,{'Amplitude';'(mV)'},'FontSize',10)
title(wave_ax,TITLE,'Interpreter','none','FontSize',10)
xlabel(spec_ax,'Time (s)','FontSize',10)