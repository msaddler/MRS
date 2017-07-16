% Script for plotting DS waveforms, spectrograms, and accel data
% MS 2017.07.09
%
% MODIFIED TO INCLUDE AUDIO-ACCEL CROSS-CORRELATION PLOTS

% Load-in savedSignals
load('/Users/mark/Documents/WHOI_SSF_2016/Mark/saddler_report/figure_data/bm15_054a_indexes_66_67_68_69_70_71_72_73_74.mat')

%call_idx = 1; % 1 and 9 (7?)
for call_idx = [9]

% cd('G:\Mark\code')
% addpath('G:\Mark\code\tools')
cd('/Users/mark/Documents/WHOI_SSF_2016/Mark/code')
addpath('/Users/mark/Documents/WHOI_SSF_2016/Mark/code/tools')

pt = 256;
pt_acc = 256;
ovl = 0.99;

FLIM = [5 250];
YLIM = [-.03 .03];

FLIM_acc = [20 120];
YLIM_acc = [-.1 .1];

% Build custom color map
MAP = 'bone';
MAP = [1:-.01:0.50, 0.49:-0.03:0 ]';
MAP = [MAP,MAP,MAP];

% Read in audio signal
savedSignal = savedSignals(call_idx,:);
fs = savedSignal{1,1}.fs;
y = savedSignal{1,1}.y;
y = y - mean(y);
t = savedSignal{1,1}.T;

% Read in accel signal
afs = savedSignal{1,2}.fs;
a1 = savedSignal{1,2}.y(:,1);
a2 = savedSignal{1,2}.y(:,2);
a3 = savedSignal{1,2}.y(:,3);
a1 = 9.81*(a1 - mean(a1));
a2 = 9.81*(a2 - mean(a2));
a3 = 8.81*(a3 - mean(a3));
ta = savedSignal{1,2}.T;

% Build figure
close all
f = figure('Position',[100 100 400 800]);
wave_ax = axes('Parent',f,'Position',[.35 .89 .55 .10]);
spec_ax = axes('Parent',f,'Position',[.35 .59 .55 .30]);
acc1_ax = axes('Parent',f,'Position',[.35 .48 .55 .10]);
acc2_ax = axes('Parent',f,'Position',[.35 .38 .55 .10]);
acc3_ax = axes('Parent',f,'Position',[.35 .28 .55 .10]);
xcor_ax = axes('Parent',f,'Position',[.35 .03 .55 .18]);

% Plot audio
MS_waveform(wave_ax,t,y,YLIM);
MS_spectrogram(spec_ax,y,fs,pt,ovl,FLIM,MAP);

% Plot accel channels (waveforms?)
MS_waveform(acc1_ax,ta,a1,YLIM_acc);
MS_waveform(acc2_ax,ta,a2,YLIM_acc);
MS_waveform(acc3_ax,ta,a3,YLIM_acc);

% Plot accel channels (spectrograms?)
% MS_spectrogram(acc1_ax,a1,afs,pt_acc,ovl,FLIM_acc,MAP);
% MS_spectrogram(acc2_ax,a2,afs,pt_acc,ovl,FLIM_acc,MAP);
% MS_spectrogram(acc3_ax,a3,afs,pt_acc,ovl,FLIM_acc,MAP);

% Link all time axes
%linkaxes([spec_ax, wave_ax, acc1_ax, acc2_ax, acc3_ax],'x')
% Link colormap axes (if applicable)
CLIM = get(acc1_ax,'clim');
set(acc2_ax,'clim',CLIM)
set(acc3_ax,'clim',CLIM)
% Clear uneccessary x-axis labels
set([acc1_ax, acc2_ax, wave_ax, spec_ax],'xticklabel',[])
set(spec_ax,'ytick',[0:50:200])
set(wave_ax,'ytick',[-.020 0 .020])
set(acc1_ax,'ytick',[-.05 0 .05])
set(acc2_ax,'ytick',[-.05 0 .05])
set(acc3_ax,'ytick',[-.05 0 .05])

ylabel(spec_ax,{'Frequency';'(Hz)'})
ylabel(wave_ax,{'Audio';'(mV)'})
ylabel(acc1_ax,{'';'X'})
ylabel(acc2_ax,{'Acceleration (ms^{-2})';'Y'})
ylabel(acc3_ax,{'';'Z'})
xlabel(acc3_ax,'Time (s)')

set([spec_ax,wave_ax,acc1_ax,acc2_ax,acc3_ax, xcor_ax],'FontSize',12)

% Get signals for audio-accel cross-correlation
O_dec = savedSignal{1,1}.y_dec;
A_dec = savedSignal{1,1}.A_dec;
xcor_list = cell(1,3);
lags_list = cell(1,3);

% Calculate audio-accel cross-correlations
for accIdx = 1:3
    A_sig = A_dec(:, accIdx);
    O_sig = O_dec;
    
    % Normalize signal 1
    A_sig = (A_sig - mean(A_sig))/std(A_sig);
    % Normalize signal 2
    O_sig = (O_sig - mean(O_sig))/std(O_sig);
    
    % Calculate cross-correlation, take abs(), normalize by lengths
    [tmp_xcor, tmp_lag] = xcorr(A_sig,O_sig);
    tmp_xcor = abs(tmp_xcor);
    tmp_xcor = tmp_xcor / (sqrt(length(A_sig))*sqrt(length(O_sig)));
    
    xcor_list{accIdx} = tmp_xcor;
    lags_list{accIdx} = tmp_lag;
end

hold(xcor_ax, 'on')
xh = plot(xcor_ax, lags_list{1} .* 1/afs, xcor_list{1}, 'm');
yh = plot(xcor_ax, lags_list{2} .* 1/afs, xcor_list{2}, 'g');
zh = plot(xcor_ax, lags_list{3} .* 1/afs, xcor_list{3}, 'b');
set(xcor_ax, 'ylim', [0, 1])
ylabel(xcor_ax, {'Cross-Correlation';'(Normalized)'})

xstr = ['X, coef = ', num2str(max(xcor_list{1}), 2)];
ystr = ['Y, coef = ', num2str(max(xcor_list{2}), 2)];
zstr = ['Z, coef = ', num2str(max(xcor_list{3}), 2)];

legend(xstr, ystr, zstr, 'location', 'northeast')

disp(call_idx)
disp('.')

end