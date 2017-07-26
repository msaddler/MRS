function signalParams = dtagcorrelate(savedSignals,varargin)
% MS 2016.08.10
%
% Function to analyze low-frequency baleen whale calls and use
% cross-correlation methods for phase analysis of accelerometer and
% acoustic signals (ideas from Dan Zitterbart). Variable input argument
% must be 'save' inorder for image files to be saved.
%
% Inputs:
%   savedSignals: cell array of signal structures created by dtagsignal.m
%                 NOTE: signal structures must include decimated audio and
%                 accelerometer time series (y_dec and A_dec); see
%                 decimate_for_correlation.m
%   varargin: 'save' in order to save all image files
%
% Outputs:
%   signalParams: structure with fields containing cell arrays for each
%                 set of signal parameter values
%

% Define directory for saving figures
saveDir = '/Users/mark/Documents/WHOI_SSF_2016/Mark/data_bm/dtag_analysis/dtaganalyze_images/';
% Define filename for temporary output file containing signal parameters
tmpSaveStr = '/Users/mark/Documents/WHOI_SSF_2016/Mark/data_bm/TEMP.mat';

% Initialize output variables
signalParams = struct();

% Fourier Analysis Parameters
pt_aud = 128; % FFT size (audio data)
pt_acc = 128; % FFT size (accelerometer data)
ovl = 0.99; % FFT overlap
FLIM_AUDIO = [0 250]; % Limits for audio spectrogram
MAP = 'bone'; % Spectrogam color scheme

% Prepare figure and axes to save as image files
close all
f = figure('Position',[1950 -50 1900 1000]);
fpan = uipanel('Parent',f,'Position',[.005 .005 .990 .990]);
ax = gobjects(4,7);
for row = 1:4
    for col = 1:7
        w = 0.12;
        h = 0.20;
        x = 0.05 + (col-1)*0.135;
        y = 1 - 0.24*row;
        ax(row,col) = axes('Parent',fpan,'Position',[x y w h]);
    end
end

% Iterate through all saved signals
for i = 1:size(savedSignals,1)
    
    %%% ACOUSTIC AND ACCELEROMETER SIGNAL CROSS-CORRELATIONS %%%
    O_dec = savedSignals{i,1}.y_dec;
    A_dec = savedSignals{i,1}.A_dec;
    
    % sigs = {audio,acc1,acc2,acc3}
    sigs{1} = O_dec;
    sigs{2} = A_dec(:,1);
    sigs{3} = A_dec(:,2);
    sigs{4} = A_dec(:,3);
    
    % Calculate signal cross-correlations, coefficients, and lags
    sigs_xcor = cell(4,4);
    sigs_coef = cell(4,4);
    sigs_lags = cell(4,4);
    for row = 1:4
        for col = 1:4
            tmp_sigRow = sigs{row};
            tmp_sigCol = sigs{col};
            
            % Normalize signal 1
            tmp_sigRow = (tmp_sigRow - mean(tmp_sigRow))/std(tmp_sigRow);
            % Normalize signal 2
            tmp_sigCol = (tmp_sigCol - mean(tmp_sigCol))/std(tmp_sigCol);
            % Calculate cross-correlation, take abs(), normalize by lengths
            [tmp_xcor, tmp_lag] = xcorr(tmp_sigRow,tmp_sigCol);
            tmp_xcor = abs(tmp_xcor);
            tmp_xcor = tmp_xcor / (sqrt(length(tmp_sigRow))*sqrt(length(tmp_sigCol)));
            
            [sigs_coef{row,col}, tmp_idx] = max(tmp_xcor);
            sigs_lags{row,col} = tmp_lag(tmp_idx);
            sigs_xcor{row,col} = tmp_xcor;
            
            % Plot cross-correlations
            plot(ax(row,col+3),tmp_lag,tmp_xcor)
        end
    end
    
    % Format and annotate the cross-correlation axes
    set(ax(1:4,4:7),'xtick',[-500 0 500],'color','none','ylim',[0 1])
    linkaxes(ax(1:4,4:7),'xy')
    for row = 1:4
        for col = 1:4
            MS_annotate_axes(ax(row,col+3), ...
                ['coef =',num2str(sigs_coef{row,col}), ...
                '  delay = ',num2str(sigs_lags{row,col})])
        end
    end
    
    %%% ACOUSTIC SIGNAL PARAMETER CALCULATIONS %%%
    
    % Get audio signal and audio noise
    y = savedSignals{i,1}.y;
    y_noise = [savedSignals{i,1}.noise_pre;savedSignals{i,1}.noise_post];
    % Subtract by means to adjust for zero-offset
    y = y-mean(y,1);
    y_noise = y_noise-mean(y_noise,1);
    % Calculate audio SNR, RMS, Peak-to-Peak, and 97% energy duration
    [aud_SNR,aud_RMS,aud_P2P,aud_dur] = ...
        amplitude_calculations(y,y_noise);
    aud_dur = aud_dur/savedSignals{i,1}.fs; % Convert dur to seconds
    % Calculate frequency parameters for audio signal
    [aud_fpeak,aud_fcenter,aud_BW3dB,aud_BW10dB] = ...
        frequency_calculations(y,y_noise,savedSignals{i,1}.fs,pt_aud,ovl);
    
    %%% ACCEL SIGNAL PARAMETER CALCULATIONS %%%
    
    % Get accelerometer signal and accelerometer noise
    A = savedSignals{i,2}.y;
    A_noise = [savedSignals{i,2}.noise_pre;savedSignals{i,2}.noise_post];
    % Subtract by means to adjust for zero-offset
    for itrA = 1:size(A,2) % Iterate across all accelerometer channels
        A(:,itrA) = A(:,itrA) - mean(A(:,itrA));
        A_noise(:,itrA) = A_noise(:,itrA) - mean(A_noise(:,itrA));
    end
    % Calculate accel SNR, RMS, Peak-to-Peak, and 97% energy duration
    [ac1_SNR,ac1_RMS,ac1_P2P,ac1_dur] = ...
        amplitude_calculations(A(:,1),A_noise(:,1)); % axis-1
    [ac2_SNR,ac2_RMS,ac2_P2P,ac2_dur] = ...
        amplitude_calculations(A(:,2),A_noise(:,2)); % axis-2
    [ac3_SNR,ac3_RMS,ac3_P2P,ac3_dur] = ...
        amplitude_calculations(A(:,3),A_noise(:,3)); % axis-3
    % Convert calculated accel signal durations to seconds
    ac1_dur = ac1_dur/savedSignals{i,2}.fs;
    ac2_dur = ac2_dur/savedSignals{i,2}.fs;
    ac3_dur = ac3_dur/savedSignals{i,2}.fs;
    % Calculate frequency parameters for accel signal
    [ac1_fpeak,ac1_fcenter,ac1_BW3dB,ac1_BW10dB] = ...
        frequency_calculations(A(:,1),A_noise(:,1),savedSignals{i,2}.fs,pt_acc,ovl); % axis-1
    [ac2_fpeak,ac2_fcenter,ac2_BW3dB,ac2_BW10dB] = ...
        frequency_calculations(A(:,2),A_noise(:,2),savedSignals{i,2}.fs,pt_acc,ovl); % axis-2
    [ac3_fpeak,ac3_fcenter,ac3_BW3dB,ac3_BW10dB] = ...
        frequency_calculations(A(:,3),A_noise(:,3),savedSignals{i,2}.fs,pt_acc,ovl); % axis-3
    
    %%% STORE PARAMETERS IN OUTPUT DATA STRUCTURES %%%
    
    % Signal-identifying information
    signalParams.file{i,1} = savedSignals{i,1}.file; % FILENAME
    signalParams.datestr{i,1} = datestr(savedSignals{i,1}.T(1),'yyyy.mm.dd HH:MM:SS'); % DATESTRING
    signalParams.call{i,1} = savedSignals{i,1}.call; % CALL TYPE
    signalParams.comment{i,1} = savedSignals{i,1}.comment; % COMMENT
    signalParams.time{i,1} = savedSignals{i,1}.T(1); % CALL TIME (matlab datenum)
    signalParams.dur{i,1} = savedSignals{i,1}.T(end)-savedSignals{i,1}.T(1); % CALL DURATION (matlab datenum)
    signalParams.duration97{i,1} = aud_dur; % 97% energy duration
    signalParams.depth{i,1} = mean(savedSignals{i,2}.p); % Depth (m)
    
    % Signal frequency measurements
    signalParams.aud_peakFreq{i,1} = aud_fpeak;
    signalParams.aud_centFreq{i,1} = aud_fcenter;
    signalParams.aud_startFreq{i,1} = NaN; % To be calculated from contour
    signalParams.aud_endFreq{i,1} = NaN; % To be calculated from contour
    signalParams.aud_maxFreq{i,1} = NaN; % To be calculated from contour
    signalParams.aud_minFreq{i,1} = NaN; % To be calculated from contour
    % Signal amplitude measurements
    signalParams.aud_SNR{i,1} = aud_SNR;
    signalParams.ac1_SNR{i,1} = ac1_SNR;
    signalParams.ac2_SNR{i,1} = ac2_SNR;
    signalParams.ac3_SNR{i,1} = ac3_SNR;
    signalParams.aud_RMS{i,1} = aud_RMS;
    signalParams.ac1_RMS{i,1} = ac1_RMS;
    signalParams.ac2_RMS{i,1} = ac2_RMS;
    signalParams.ac3_RMS{i,1} = ac3_RMS;
    signalParams.aud_P2P{i,1} = aud_P2P;
    signalParams.ac1_P2P{i,1} = ac1_P2P;
    signalParams.ac2_P2P{i,1} = ac2_P2P;
    signalParams.ac3_P2P{i,1} = ac3_P2P;
    % Signal phase measurements
    signalParams.ox_coef{i,1} = sigs_coef{1,2};
    signalParams.oy_coef{i,1} = sigs_coef{1,3};
    signalParams.oz_coef{i,1} = sigs_coef{1,4};
    signalParams.ox_delay{i,1} = sigs_lags{1,2};
    signalParams.oy_delay{i,1} = sigs_lags{1,3};
    signalParams.oz_delay{i,1} = sigs_lags{1,4};
    signalParams.xy_coef{i,1} = sigs_coef{2,3};
    signalParams.xz_coef{i,1} = sigs_coef{2,4};
    signalParams.yz_coef{i,1} = sigs_coef{3,4};
    signalParams.xy_delay{i,1} = sigs_lags{2,3};
    signalParams.xz_delay{i,1} = sigs_lags{2,4};
    signalParams.yz_delay{i,1} = sigs_lags{3,4};
    signalParams.oo_coef{i,1} = sigs_coef{1,1};
    signalParams.xx_coef{i,1} = sigs_coef{2,2};
    signalParams.yy_coef{i,1} = sigs_coef{3,3};
    signalParams.zz_coef{i,1} = sigs_coef{4,4};
    signalParams.oo_delay{i,1} = sigs_lags{1,1};
    signalParams.xx_delay{i,1} = sigs_lags{2,2};
    signalParams.yy_delay{i,1} = sigs_lags{3,3};
    signalParams.zz_delay{i,1} = sigs_lags{4,4};
    
    %%% GENERATE FIGURE TO SAVE AS IMAGE FILE %%%
    
    % Generate plots to save as .png file
    titleStr = ['File: ',savedSignals{i,1}.file,...
        '      Time: ',datestr(savedSignals{i,1}.T(1),'yyyy.mm.dd HH:MM:SS'),...
        '      Call type: ',savedSignals{i,1}.call];
    set(fpan,'Title',titleStr)
    
    % Plot all signal waveforms
    MS_waveform(ax(1,1),[1:length(O_dec)]/savedSignals{i,1}.fs_dec,O_dec); %MS_waveform(ax(1,1),savedSignals{i,1}.T,y);
    MS_waveform(ax(2,1),[1:length(A_dec)]/savedSignals{i,1}.fs_dec,A_dec(:,1)); %MS_waveform(ax(2,1),savedSignals{i,2}.T,A(:,1));
    MS_waveform(ax(3,1),[1:length(A_dec)]/savedSignals{i,1}.fs_dec,A_dec(:,2)); %MS_waveform(ax(3,1),savedSignals{i,2}.T,A(:,2));
    MS_waveform(ax(4,1),[1:length(A_dec)]/savedSignals{i,1}.fs_dec,A_dec(:,3)); %MS_waveform(ax(4,1),savedSignals{i,2}.T,A(:,3));
    
    % Plot all signal spectrograms
    [S_aud,F_aud,T_aud] = MS_spectrogram(ax(1,2),O_dec,savedSignals{i,1}.fs_dec,pt_aud,ovl,savedSignals{i,2}.bandpass,MAP);
    MS_spectrogram(ax(2,2),A_dec(:,1),savedSignals{i,2}.fs,pt_acc,ovl,savedSignals{i,2}.bandpass,MAP);
    MS_spectrogram(ax(3,2),A_dec(:,2),savedSignals{i,2}.fs,pt_acc,ovl,savedSignals{i,2}.bandpass,MAP);
    MS_spectrogram(ax(4,2),A_dec(:,3),savedSignals{i,2}.fs,pt_acc,ovl,savedSignals{i,2}.bandpass,MAP);
    
    % Plot all signal power spectra
    noiseIdx = 0.15*savedSignals{i,1}.fs_dec;
    MS_powerspectrum_noise_correction(ax(1,3),O_dec,O_dec([1:noiseIdx,end-noiseIdx:end]),savedSignals{i,1}.fs_dec,pt_aud,ovl,FLIM_AUDIO,MAP);
    MS_powerspectrum_noise_correction(ax(2,3),A_dec(:,1),A_dec([1:noiseIdx,end-noiseIdx:end],1),savedSignals{i,2}.fs,pt_acc,ovl,FLIM_AUDIO);
    MS_powerspectrum_noise_correction(ax(3,3),A_dec(:,2),A_dec([1:noiseIdx,end-noiseIdx:end],2),savedSignals{i,2}.fs,pt_acc,ovl,FLIM_AUDIO);
    MS_powerspectrum_noise_correction(ax(4,3),A_dec(:,3),A_dec([1:noiseIdx,end-noiseIdx:end],3),savedSignals{i,2}.fs,pt_acc,ovl,FLIM_AUDIO);
    
    % Format Plots
    title(ax(1,1),'Waveform','Fontsize',8,'Fontweight','normal')
    title(ax(1,2),'Spectrogram','Fontsize',8,'Fontweight','normal')
    title(ax(1,3),'Power Spectrum','Fontsize',8,'Fontweight','normal')
    ylabel(ax(1,1),'Audio','Fontsize',8,'Fontweight','normal')
    ylabel(ax(2,1),'Accel (X)','Fontsize',8,'Fontweight','normal')
    ylabel(ax(3,1),'Accel (Y)','Fontsize',8,'Fontweight','normal')
    ylabel(ax(4,1),'Accel (Z)','Fontsize',8,'Fontweight','normal')
    title(ax(1,4),'Audio','Fontsize',8,'Fontweight','normal')
    title(ax(1,5),'Accel (X)','Fontsize',8,'Fontweight','normal')
    title(ax(1,6),'Accel (Y)','Fontsize',8,'Fontweight','normal')
    title(ax(1,7),'Accel (Z)','Fontsize',8,'Fontweight','normal')
    linkaxes(ax(2:4,1),'xy')
    
    % Format and annotate axes with call parameters
    linkaxes(ax(:,1),'xy')
    set(ax(:,1),'ylim',[-0.04 0.04])
    MS_annotate_axes(ax(1,1),{['SNR = ',num2str(aud_SNR)];['peak-to-peak = ',num2str(aud_P2P)]});
    MS_annotate_axes(ax(2,1),{['SNR = ',num2str(ac1_SNR)];['peak-to-peak = ',num2str(ac1_P2P)]});
    MS_annotate_axes(ax(3,1),{['SNR = ',num2str(ac2_SNR)];['peak-to-peak = ',num2str(ac2_P2P)]});
    MS_annotate_axes(ax(4,1),{['SNR = ',num2str(ac3_SNR)];['peak-to-peak = ',num2str(ac3_P2P)]});
    
    %%% SAVE FIGURE AS IMAGE FILE AND SAVE BACKUP FILES %%%
    
    if ~isempty(varargin) && strcmp(varargin{1},'save')
        % If variable input argument is 'save', then save each figure as
        % .png image file and save temporary backup signalParams file.
        
        % Save plots as .png image
        set(f,'PaperPositionMode','auto')
        %         saveStr = [saveDir,'accZcoef_',num2str(round(sigs_coef{1,4},2)),'_',...
        %             savedSignals{i,1}.call,'_',datestr(savedSignals{i,1}.T(1),'yyyy-mm-dd_HH-MM-SS'),'.png'];
        saveStr = [saveDir,'index',num2str(i,'%02.0f'),'_',savedSignals{i,1}.call,'_',datestr(savedSignals{i,1}.T(1),'yyyy-mm-dd_HH-MM-SS'),'.png'];
        print(f,saveStr,'-dpng','-r300')
        
        % Save output structures at each step to minimize lost time
        save(tmpSaveStr,'signalParams');
    end

end

disp('COMPLETE!')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS CALLED BY DRIVER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [SNR,RMS,P2P,dur] = amplitude_calculations(signal,noise,varargin)
        % Function to calculate signal amplitude parameters
        % (modified from code provided by Max Kaplan)
        % Variable input argument allows user to specific energy bounds of
        % signal to analyze (default is middle 97% of signal energy)
        
        if isempty(varargin)
            E_trim = 0.97; % Default middle 97% of energy
        else
            E_trim = varargin{1,1}; % User input energy fraction
        end
        trim = (1-E_trim)/2;
        
        % Calculate mean-square amplitude of noise
        ms_noise = (1/length(noise))*sum(noise.^2);
        
        % Trim signal at energy bounds (default bounds 97% energy)
        E_rel = sum(signal.^2); % relative energy (sum of squared pressures)
        E_cum = cumsum(signal.^2)-ms_noise; % cumulative sum of energy (cumulative sum of squared pressures)
        E_start = find(E_cum > trim*E_rel,1);
        E_end = find(E_cum < (1-trim)*E_rel,1,'last');
        signal = signal(E_start:E_end);
        dur = length(signal); % Return the duration of the energy-bounded signal in # of samples
        
        % Calculate mean-square amplitude of signal (bounded by 97% energy)
        ms_signal = (1/length(signal))*sum(signal.^2);
        % Calculate signal-to-noise ratio in dB
        if ms_signal-ms_noise >= 0
            SNR = 10*log10((ms_signal-ms_noise)/ms_noise);
        else
            SNR = NaN;
        end
        
        % Calculate signal root mean-square and peak-to-peak amplitudes
        RMS = rms(signal);
        P2P = peak2peak(signal);
    end

    function [f_peak,f_center,BW_3dB,BW_10dB]...
            = frequency_calculations(signal,noise,fs,pt,ovl,varargin)
        % Function to calculate signal frequency parameters
        % (modified from code provided by Max Kaplan)
        % Variable input argument allows user to specific energy bounds of
        % signal to analyze (default is middle 97% of signal energy)
        
        if isempty(varargin)
            E_trim = 0.97; % Default middle 97% of energy
        else
            E_trim = varargin{1,1}; % User input energy fraction
        end
        trim = (1-E_trim)/2;
        
        % Calculate mean-square amplitude of noise
        ms_noise = (1/length(noise))*sum(noise.^2);
        
        % Trim signal at energy bounds (default bounds 97% energy)
        E_rel = sum(signal.^2); % relative energy (sum of squared pressures)
        E_cum = cumsum(signal.^2)-ms_noise; % cumulative sum of energy (cumulative sum of squared pressures)
        E_start = find(E_cum > trim*E_rel,1);
        E_end = find(E_cum < (1-trim)*E_rel,1,'last');
        signal = signal(E_start:E_end); % Trim signal to energy bounds
        
        % Set up Fourier analysis parameters
        window = hamming(pt);
        noverlap = floor(ovl*pt);
        if length(signal) < pt % Work-around to avoid short segments causing errors
            warning('Short segment: Fourier parameters modified')
            window = hamming(length(signal));
            noverlap =  floor(ovl*length(signal));
        end
        
        % Calculate peak and center frequencies from periodogram
        [Pxx,Fxx]=pwelch(signal,window,noverlap,pt,fs,'onesided');
        [~,peakIdx] = max(Pxx); % Find where power spectrum has maximum
        f_peak = Fxx(peakIdx); % Peak frequency
        f_center = sum(Pxx.*Fxx)/sum(Pxx); % Center frequency
        
        % Calculate -3dB and -10db bandwidths
        PxxdB = 10*log10(Pxx); % Convert power spectrum to dB
        ampPeak = PxxdB(peakIdx); % Peak amplitude in dB
        
        BW = -3; % Calculate -3dB bandwidth
        idx = find(PxxdB-ampPeak > BW);
        BW_3dB = [min(Fxx(idx)), max(Fxx(idx))];
        
        BW = -10; % Calculate -10dB bandwidth
        idx = find(PxxdB-ampPeak > BW);
        BW_10dB = [min(Fxx(idx)), max(Fxx(idx))];
    end

end