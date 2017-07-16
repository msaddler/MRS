% dtaganalyze_bandpower.m
% MS 2016.07.29

% Load-in the following variables:
%   signalParams
%   signalContours
%   savedSignals

% PARAMETERS
pt_aud = 256; % FFT size (audio)
pt_acc = 128; % FFT size (accel)
ovl = 0.99; % Fraction of overlap for FFT
FLIM_AUDIO = [0 500]; % Spectrogram frequency limits (audio)
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
close all
f = figure('Position',[2050 -50 1600 1000]);
fpan = uipanel('Parent',f,'Position',[.01 .01 .98 .98]);
ax = gobjects(4,3);
for row = 1:4
    for col = 1:3
        w = 0.27;
        h = 0.20;
        x = 0.06 + (col-1)*0.30;
        y = 1 - 0.24*row;
        ax(row,col) = axes('Parent',fpan,'Position',[x y w h]);
    end
end
linkaxes(ax(2:4,1),'xy')
linkaxes(ax(1:4,3),'xy')

% Iterate through all saved signals
for i = 1:size(savedSignals,1)
    % Get audio signal and audio noise
    y = savedSignals{i,1}.y;
    yn = [savedSignals{i,1}.noise_pre; savedSignals{i,1}.noise_post];
    % Subtract by means to adjust for zero-offset
    y = y-mean(y,1);
    yn = yn-mean(yn,1);
    % Calculate signal-to-noise ratio
    snr_y = MS_snr(y,yn);
    
    % Get accelerometer signal and accelerometer noise
    A = savedSignals{i,2}.y;
    An = [savedSignals{i,2}.noise_pre; savedSignals{i,2}.noise_post];
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
    
    % Display identifying information for signal on the plot panel
    titleStr = ['File: ',savedSignals{i,1}.file,...
        '      Time: ',datestr(savedSignals{i,1}.T(1),'yyyy.mm.dd HH:MM:SS'),...
        '      Call type: ',savedSignals{i,1}.call];
    set(fpan,'Title',titleStr)
    
    % Plot all signal waveforms
    MS_waveform(ax(1,1),savedSignals{i,1}.T,y);
    MS_waveform(ax(2,1),savedSignals{i,2}.T,A(:,1));
    MS_waveform(ax(3,1),savedSignals{i,2}.T,A(:,2));
    MS_waveform(ax(4,1),savedSignals{i,2}.T,A(:,3));
    
    % Plot all signal spectrograms
    MS_spectrogram(ax(1,2),y,savedSignals{i,1}.fs,pt_aud,ovl,FLIM_AUDIO,MAP);
    MS_spectrogram(ax(2,2),A(:,1),savedSignals{i,2}.fs,pt_acc,ovl,savedSignals{i,2}.bandpass,MAP);
    MS_spectrogram(ax(3,2),A(:,2),savedSignals{i,2}.fs,pt_acc,ovl,savedSignals{i,2}.bandpass,MAP);
    MS_spectrogram(ax(4,2),A(:,3),savedSignals{i,2}.fs,pt_acc,ovl,savedSignals{i,2}.bandpass,MAP);
    
    % Plot all signal power spectra
    [Pxx,Fxx,Pyy,Fyy] = MS_powerspectrum_noise_correction(ax(1,3),y,yn,savedSignals{i,1}.fs,pt_aud,ovl,FLIM_AUDIO,MAP);
    [Pxx1,Fxx1,Pyy1,Fyy1] = MS_powerspectrum_noise_correction(ax(2,3),A(:,1),An(:,1),savedSignals{i,2}.fs,pt_acc,ovl,FLIM_AUDIO);
    [Pxx2,Fxx2,Pyy2,Fyy2] = MS_powerspectrum_noise_correction(ax(3,3),A(:,2),An(:,2),savedSignals{i,2}.fs,pt_acc,ovl,FLIM_AUDIO);
    [Pxx3,Fxx3,Pyy3,Fyy3] = MS_powerspectrum_noise_correction(ax(4,3),A(:,3),An(:,3),savedSignals{i,2}.fs,pt_acc,ovl,FLIM_AUDIO);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % POWER SPECTRA CALCULATIONS (BELOW)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Determine the audio signal's frequency band from signalContours
    f_band(1) = min(signalContours.f_audio{i});
    f_band(2) = max(signalContours.f_audio{i});
    
    % Convert all audio and accel signal and noise spectra to dB
    Pxx = 10*log10(abs(Pxx)); Pyy = 10*log10(abs(Pyy));
    Pxx1 = 10*log10(abs(Pxx1)); Pyy1 = 10*log10(abs(Pyy1));
    Pxx2 = 10*log10(abs(Pxx2)); Pyy2 = 10*log10(abs(Pyy2));
    Pxx3 = 10*log10(abs(Pxx3)); Pyy3 = 10*log10(abs(Pyy3));
    
    % Calculate peak frequency
    [~,idx] = max(Pxx);
    f_peak = Fxx(idx);
    
    % Calculate audio signal's power in frequency band
    idx = find(Fxx >= f_band(1) & Fxx <= f_band(2));
    audio_bandpower = Pxx(idx)-Pyy(idx);
    audio_bandpower = sum(audio_bandpower);%(find(audio_bandpower>0)));
    % Calculate accel 1 signal's power in frequency band
    idx = find(Fxx1 >= f_band(1) & Fxx1 <= f_band(2));
    acc1_bandpower = Pxx1(idx)-Pyy1(idx);
    acc1_bandpower = sum(acc1_bandpower);%(find(acc1_bandpower>0)));
    % Calculate accel 2 signal's power in frequency band
    idx = find(Fxx2 >= f_band(1) & Fxx2 <= f_band(2));
    acc2_bandpower = Pxx2(idx)-Pyy2(idx);
    acc2_bandpower = sum(acc2_bandpower);%(find(acc2_bandpower>0)));
    % Calculate accel 3 signal's power in frequency band
    idx = find(Fxx3 >= f_band(1) & Fxx3 <= f_band(2));
    acc3_bandpower = Pxx3(idx)-Pyy3(idx);
    acc3_bandpower = sum(acc3_bandpower);%(find(acc3_bandpower>0)));
    
    % Calculate ratios of accelerometer to audio bandpower
    ratio1 = acc1_bandpower/audio_bandpower;
    ratio2 = acc2_bandpower/audio_bandpower;
    ratio3 = acc3_bandpower/audio_bandpower;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % POWER SPECTRA CALCULATIONS (ABOVE)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Format plots
    title(ax(1,1),'Waveform','Fontsize',8,'Fontweight','normal')
    title(ax(1,2),'Spectrogram','Fontsize',8,'Fontweight','normal')
    title(ax(1,3),'Power Spectrum','Fontsize',8,'Fontweight','normal')
    ylabel(ax(1,1),'Audio','Fontsize',8,'Fontweight','normal')
    ylabel(ax(2,1),'Accel (X)','Fontsize',8,'Fontweight','normal')
    ylabel(ax(3,1),'Accel (Y)','Fontsize',8,'Fontweight','normal')
    ylabel(ax(4,1),'Accel (Z)','Fontsize',8,'Fontweight','normal')
    
    % Display parameters on plots
    MS_annotate_axes(ax(1,1),['SNR=',num2str(snr_y)]);
    MS_annotate_axes(ax(2,1),['SNR=',num2str(snr_A(1))]);
    MS_annotate_axes(ax(3,1),['SNR=',num2str(snr_A(2))]);
    MS_annotate_axes(ax(4,1),['SNR=',num2str(snr_A(3))]);
    
    MS_annotate_axes(ax(1,3),{['Peak freq =',num2str(f_peak)];...
        ['Start freq =',num2str(signalContours.f_audio{i}(1))];...
        ['End freq =',num2str(signalContours.f_audio{i}(end))];...
        ['Bandpower =',num2str(audio_bandpower)]});
    MS_annotate_axes(ax(2,3),{['Bandpower =',num2str(acc1_bandpower)];...
        ['Bandpower ratio =',num2str(ratio1)];...
        [num2str(f_band(1)),'Hz to ',num2str(f_band(2)),'Hz']});
    MS_annotate_axes(ax(3,3),{['Bandpower =',num2str(acc2_bandpower)];...
        ['Bandpower ratio =',num2str(ratio2)];...
        [num2str(f_band(1)),'Hz to ',num2str(f_band(2)),'Hz']});
    MS_annotate_axes(ax(4,3),{['Bandpower =',num2str(acc3_bandpower)];...
        ['Bandpower ratio =',num2str(ratio3)];...
        [num2str(f_band(1)),'Hz to ',num2str(f_band(2)),'Hz']});
    
    % Add parameters to 'signalParams' output structure
    signalParams.audio_bandpower{i,1} = audio_bandpower;
    signalParams.acc1_bandpower{i,1} = acc1_bandpower;
    signalParams.acc2_bandpower{i,1} = acc2_bandpower;
    signalParams.acc3_bandpower{i,1} = acc3_bandpower;
    
    % Save plots as .png image
    avgRatio = mean([ratio1,ratio2,ratio3]);
    nonfocalIndex = round(100*(1-avgRatio),1);
    saveStr = [saveDir,num2str(nonfocalIndex),'_',...
        savedSignals{i,1}.call,'_',datestr(savedSignals{i,1}.T(1),'yyyy-mm-dd_HH-MM-SS'),'.png'];
    set(f,'PaperPositionMode','auto')
    print(f,saveStr,'-dpng','-r300')
end
disp('COMPLETE')
clearvars -except savedSignals signalContours signalParams