% Script for plotting spectral analysis figures of downsweep calls
% MS 2016.08.02

% modified from dtaganalyze_downsweep.m

% Load-in the following variables:
%   signalParams
%   signalContours
%   savedSignals

% PARAMETERS
pt_aud = 256; % FFT size (audio)
pt_acc = 128; % FFT size (accel)
ovl = 0.99; % Fraction of overlap for FFT
FLIM_AUDIO = [0 250]; % Spectrogram frequency limits (audio)
FLIM_ACCEL = [0 250]; % Spectrogram frequency limits (accel)
YLIM_AUDIO = [0 0];
YLIM_ACCEL = [0 0];

% Build custom color map for spectrograms
MAP = [1:-.01:0.50, 0.49:-0.03:0 ]';
MAP = [MAP,MAP,MAP];

% Define directory to save images
saveDir = 'G:\Mark\data_bm\dtag_analysis\dtaganalyze_images\';

% Ensure all called functions are in the path
cd('G:\Mark\code')
addpath('G:\Mark\code\tools')

% Set up figure
%close all
f = figure('Position',[2050 -115 750 1100]);
ax = gobjects(4,2);
for row = 1:4
    for col = 1:2
        w = 0.35;
        h = 0.20;
        x = 0.10 + (col-1)*0.50;
        y = 0.97 - 0.22*row;
        ax(row,col) = axes('Parent',f,'Position',[x y w h]);
    end
end
linkaxes(ax(2:4,1),'xy')
linkaxes(ax(1:4,2),'xy')

% Get audio signal and audio noise
y = savedSignals{1,1}.y;
yn = [savedSignals{1,1}.noise_pre; savedSignals{1,1}.noise_post];
% Subtract by means to adjust for zero-offset
y = y-mean(y,1);
yn = yn-mean(yn,1);
% Calculate signal-to-noise ratio
snr_y = MS_snr(y,yn);

% Get accelerometer signal and accelerometer noise
A = savedSignals{1,2}.y;
An = [savedSignals{1,2}.noise_pre; savedSignals{1,2}.noise_post];
% Subtract by means to adjust for zero-offset
for itrA = 1:size(A,2) % Iterate across all accelerometer channels
    A(:,itrA) = A(:,itrA) - mean(A(:,itrA));
    An(:,itrA) = An(:,itrA) - mean(An(:,itrA));
end
% Calculate SNR of each accel channel individually
snr_A(1) = MS_snr(A(:,1),An(:,1));
snr_A(2) = MS_snr(A(:,2),An(:,2));
snr_A(3) = MS_snr(A(:,3),An(:,3));
[snr_A_max, maxChannel] = max(snr_A);

% Plot all signal spectrograms
MS_spectrogram(ax(1,1),y,savedSignals{1,1}.fs,pt_aud,ovl,FLIM_AUDIO,MAP);
MS_spectrogram(ax(2,1),A(:,1),savedSignals{1,2}.fs,pt_acc,ovl,savedSignals{1,2}.bandpass,MAP);
MS_spectrogram(ax(3,1),A(:,2),savedSignals{1,2}.fs,pt_acc,ovl,savedSignals{1,2}.bandpass,MAP);
MS_spectrogram(ax(4,1),A(:,3),savedSignals{1,2}.fs,pt_acc,ovl,savedSignals{1,2}.bandpass,MAP);

% Plot all signal power spectra
[Pxx,Fxx,Pyy,Fyy] = MS_powerspectrum_noise_correction(ax(1,2),y,yn,savedSignals{1,1}.fs,pt_aud,ovl,FLIM_AUDIO,MAP);
[Pxx1,Fxx1,Pyy1,Fyy1] = MS_powerspectrum_noise_correction(ax(2,2),A(:,1),An(:,1),savedSignals{1,2}.fs,pt_acc,ovl,FLIM_AUDIO);
[Pxx2,Fxx2,Pyy2,Fyy2] = MS_powerspectrum_noise_correction(ax(3,2),A(:,2),An(:,2),savedSignals{1,2}.fs,pt_acc,ovl,FLIM_AUDIO);
[Pxx3,Fxx3,Pyy3,Fyy3] = MS_powerspectrum_noise_correction(ax(4,2),A(:,3),An(:,3),savedSignals{1,2}.fs,pt_acc,ovl,FLIM_AUDIO);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POWER SPECTRA CALCULATIONS (BELOW)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the audio signal's frequency band from signalContours
f_band(1) = min(signalContours.f_audio);
f_band(2) = max(signalContours.f_audio);

% Convert all audio and accel signal and noise spectra to dB
Pxx = 10*log10(abs(Pxx)); Pyy = 10*log10(abs(Pyy));
Pxx1 = 10*log10(abs(Pxx1)); Pyy1 = 10*log10(abs(Pyy1));
Pxx2 = 10*log10(abs(Pxx2)); Pyy2 = 10*log10(abs(Pyy2));
Pxx3 = 10*log10(abs(Pxx3)); Pyy3 = 10*log10(abs(Pyy3));

% Calculate peak frequency
[~,idx] = max(Pxx);
f_peak = Fxx(idx);

% Calculate audio signal's power in frequency band (only add power if
% above the noise floor)
idx = find(Fxx >= f_band(1) & Fxx <= f_band(2));
audio_bandpower = Pxx(idx)-Pyy(idx);
audio_bandpower = sum(audio_bandpower(find(audio_bandpower>0)));
% Calculate accel 1 signal's power in frequency band (only add power if
% above the noise floor)
idx = find(Fxx1 >= f_band(1) & Fxx1 <= f_band(2));
acc1_bandpower = Pxx1(idx)-Pyy1(idx);
acc1_bandpower = sum(acc1_bandpower(find(acc1_bandpower>0)));
% Calculate accel 2 signal's power in frequency band (only add power if
% above the noise floor)
idx = find(Fxx2 >= f_band(1) & Fxx2 <= f_band(2));
acc2_bandpower = Pxx2(idx)-Pyy2(idx);
acc2_bandpower = sum(acc2_bandpower(find(acc2_bandpower>0)));
% Calculate accel 3 signal's power in frequency band (only add power if
% above the noise floor)
idx = find(Fxx3 >= f_band(1) & Fxx3 <= f_band(2));
acc3_bandpower = Pxx3(idx)-Pyy3(idx);
acc3_bandpower = sum(acc3_bandpower(find(acc3_bandpower>0)));

% Calculate ratios of accelerometer to audio bandpower
ratio1 = acc1_bandpower/audio_bandpower;
ratio2 = acc2_bandpower/audio_bandpower;
ratio3 = acc3_bandpower/audio_bandpower;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POWER SPECTRA CALCULATIONS (ABOVE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Format plots
ylabel(ax(1,1),'Audio','Fontsize',10,'Fontweight','bold')
ylabel(ax(2,1),'Accel (X)','Fontsize',10,'Fontweight','bold')
ylabel(ax(3,1),'Accel (Y)','Fontsize',10,'Fontweight','bold')
ylabel(ax(4,1),'Accel (Z)','Fontsize',10,'Fontweight','bold')
ylabel(ax(1,2),'(dB re 1 mV)','Fontsize',10,'Fontweight','normal')
ylabel(ax(2,2),'(dB re 1 mV)','Fontsize',10,'Fontweight','normal')
ylabel(ax(3,2),'(dB re 1 mV)','Fontsize',10,'Fontweight','normal')
ylabel(ax(4,2),'(dB re 1 mV)','Fontsize',10,'Fontweight','normal')

title(ax(1,1),'Spectrogram','Fontsize',10,'Fontweight','bold')
title(ax(1,2),'Power Spectrum','Fontsize',10,'Fontweight','bold')
xlabel(ax(4,1),'Time (s)','Fontsize',10,'Fontweight','normal')
xlabel(ax(4,2),'Frequency (Hz)','Fontsize',10,'Fontweight','normal')

set(ax(:,2),'ylim',[-120 -50])

set(ax(1:3,:),'xticklabel',[])
set(legend(ax(1,2)),'Visible','on')
set(legend(ax(2,2)),'Visible','off')
set(legend(ax(3,2)),'Visible','off')
set(legend(ax(4,2)),'Visible','off')

set(legend(ax(1,2)),'Location','east')
legPos = get(legend(ax(1,2)),'Position');
legPos([1,2]) = legPos([1,2]) + .02;
legPos(2) = legPos(2) + .03;
set(legend(ax(1,2)),'Position',legPos)

% Display SNR on spectrograms
MS_annotate_axes(ax(1,1),['SNR=',num2str(snr_y)]);
MS_annotate_axes(ax(2,1),['SNR=',num2str(snr_A(1))]);
MS_annotate_axes(ax(3,1),['SNR=',num2str(snr_A(2))]);
MS_annotate_axes(ax(4,1),['SNR=',num2str(snr_A(3))]);
% Display bandlimited power on spectra
MS_annotate_axes(ax(1,2),['Bandlimited power =',num2str(audio_bandpower)]);
MS_annotate_axes(ax(2,2),['Bandlimited power =',num2str(acc1_bandpower)]);
MS_annotate_axes(ax(3,2),['Bandlimited power =',num2str(acc2_bandpower)]);
MS_annotate_axes(ax(4,2),['Bandlimited power =',num2str(acc3_bandpower)]);