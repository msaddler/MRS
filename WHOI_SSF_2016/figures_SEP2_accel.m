% Script for plotting SEP2 call waveforms, spectrogram, and acceleration
% MS 2016.08.17

% Load-in savedSignals

%cd('G:\Mark\code')
%addpath('G:\Mark\code\tools')
cd('/Users/mark/Documents/WHOI_SSF_2016/Mark/code')
addpath('/Users/mark/Documents/WHOI_SSF_2016/Mark/code/tools')

pt = 256;
pt_acc = 256;
ovl = 0.99;

FLIM = [0 500];
YLIM = [-.05 .05];

FLIM_acc = [0 100];
YLIM_acc = [-.01 .01];

% Build custom color map
MAP = 'bone';
MAP = [1:-.01:0.50, 0.49:-0.03:0 ]';
MAP = [MAP,MAP,MAP];

% Read in audio signal
fs = savedSignals{1,1}.fs;
y = savedSignals{1,1}.y;
y = y - mean(y);
t = savedSignals{1,1}.T;

% Read in accel signal
afs = savedSignals{1,2}.fs;
a1 = savedSignals{1,2}.y(:,1);
a2 = savedSignals{1,2}.y(:,2);
a3 = savedSignals{1,2}.y(:,3);
a1 = a1 - mean(a1);
a2 = a2 - mean(a2);
a3 = a3 - mean(a3);
ta = savedSignals{1,2}.T;

% Build figure
close all
f = figure('Position',[50 100 1500 600]);
wave_ax = axes('Parent',f,'Position',[.10 .85 .85 .10]);
spec_ax = axes('Parent',f,'Position',[.10 .55 .85 .30]);
acc1_ax = axes('Parent',f,'Position',[.10 .42 .85 .10]);
acc2_ax = axes('Parent',f,'Position',[.10 .32 .85 .10]);
acc3_ax = axes('Parent',f,'Position',[.10 .22 .85 .10]);

% Plot audio
MS_waveform(wave_ax,t,y,YLIM);
MS_spectrogram(spec_ax,y,fs,pt,ovl,FLIM,MAP);

% Plot accel channels (waveforms?)
% MS_waveform(acc1_ax,ta,a1,YLIM_acc);
% MS_waveform(acc2_ax,ta,a2,YLIM_acc);
% MS_waveform(acc3_ax,ta,a3,YLIM_acc);

% Plot accel channels (spectrograms?)
MS_spectrogram(acc1_ax,a1,afs,pt_acc,ovl,FLIM_acc,MAP);
MS_spectrogram(acc2_ax,a2,afs,pt_acc,ovl,FLIM_acc,MAP);
MS_spectrogram(acc3_ax,a3,afs,pt_acc,ovl,FLIM_acc,MAP);

% Link all time axes
linkaxes([spec_ax, wave_ax, acc1_ax, acc2_ax, acc3_ax],'x')
% Link colormap axes (if applicable)
CLIM = get(acc1_ax,'clim');
set(acc2_ax,'clim',CLIM)
set(acc3_ax,'clim',CLIM)
% Clear uneccessary x-axis labels
set([acc1_ax, acc2_ax, wave_ax],'xticklabel',[])
set(wave_ax,'ytick',[-.03 0 .03])
set(spec_ax,'ytick',[50:100:450])
set(acc1_ax,'ytick',[25:25:75])
set(acc2_ax,'ytick',[25:25:75])
set(acc3_ax,'ytick',[25:25:75])

ylabel(spec_ax,{'Frequency';'(Hz)'})
ylabel(wave_ax,{'Amplitude';'(mV)'})
ylabel(acc1_ax,{'';'X'})
ylabel(acc2_ax,{'Accelerometer (Hz)';'Y'})
ylabel(acc3_ax,{'';'Z'})
xlabel(acc3_ax,'Time (s)')

set([spec_ax,wave_ax,acc1_ax,acc2_ax,acc3_ax],'FontSize',12)

set(spec_ax,'xticklabel',[])

% % Uncomment below to hide accelerometer axes
% delete([acc1_ax,acc2_ax,acc3_ax])
% xlabel(spec_ax,'Time (s)')
