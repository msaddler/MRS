function signalParams = dtaganalyze(savedSignals,varargin)
% MS 2016.07.14
% Script for calculating audio signal parameters from signals saved using
% dtagsignal.m

%%% PARAMETERS TO CALCULATE
% Call duration
% 95% energy duration
% 95% energy signal-to-noise ratio (SNR)
% Peak frequency, center frequency
% Start and end frequency
% -3dB and -10dB bandwidths
% RMS amplitude
% Peak-to-peak amplitude

signalParams = struct();

% Fourier Analysis Parameters
pt_aud = 256; % FFT size (audio data)
pt_acc = 64; % FFT size (accelerometer data)
ovl = 0.99; % FFT overlap
FLIM = [0 500]; % Limits for audio spectrogram

if isempty(varargin)
    saveFigs = 0;
elseif varargin{1} == 1;
    saveFigs = 1;
    saveDir = 'G:\Mark\data_bm\dtag_analysis\dtaganalyze_images\';
    
    % Prepare figure and axes to save as image files
    close all
    f = figure('Position',[2050 70 1600 900]);
    fpan = uipanel('Parent',f,'Position',[.01 .01 .98 .98]);
    ax(1,1) = axes('Parent',fpan,'Position',[.06 .04 .27 .27]);
    ax(2,1) = axes('Parent',fpan,'Position',[.06 .36 .27 .27]);
    ax(3,1) = axes('Parent',fpan,'Position',[.06 .68 .27 .27]);
    ax(1,2) = axes('Parent',fpan,'Position',[.38 .04 .27 .27]);
    ax(2,2) = axes('Parent',fpan,'Position',[.38 .36 .27 .27]);
    ax(3,2) = axes('Parent',fpan,'Position',[.38 .68 .27 .27]);
    ax(1,3) = axes('Parent',fpan,'Position',[.70 .04 .27 .27]);
    ax(2,3) = axes('Parent',fpan,'Position',[.70 .36 .27 .27]);
    ax(3,3) = axes('Parent',fpan,'Position',[.70 .68 .27 .27]);
end

% Iterate through all saved signals
for i = 1:size(savedSignals,1)
    % Get audio signal and audio noise
    y = savedSignals{i,1}.y;
    y_noise = [savedSignals{i,1}.noise_pre;savedSignals{i,1}.noise_post];
    % Subtract by means to adjust for zero-offset
    y = y-mean(y,1);
    y_noise = y_noise-mean(y_noise,1);
    % Calculate audio SNR
    [y_SNR, dur] = snr_calc(y,y_noise);
    
    % Get accelerometer signal and accelerometer noise
    A = savedSignals{i,2}.y;
    A_noise = [savedSignals{i,2}.noise_pre;savedSignals{i,2}.noise_post];
    % Subtract by means to adjust for zero-offset
    for itrA = 1:size(A,2) % Iterate across all accelerometer channels
        A(:,itrA) = A(:,itrA) - mean(A(:,itrA));
        A_noise(:,itrA) = A_noise(:,itrA) - mean(A_noise(:,itrA));
    end
    % Calculate SNR of each accel channel individually
    [A_SNR(1),~] = snr_calc(A(:,1),A_noise(:,1));
    [A_SNR(2),~] = snr_calc(A(:,2),A_noise(:,2));
    [A_SNR(3),~] = snr_calc(A(:,3),A_noise(:,3));
    [A_max_SNR, maxChannel] = max(A_SNR);
    % Calculate SNR of mean across all accel channels
    A_mean = mean(A,2);
    A_mean_noise = mean(A_noise,2);
    [A_mean_SNR,~] = snr_calc(A_mean,A_mean_noise);
    
    % Calculate frequency parameters for audio signal
    [t_start,t_end,f_start,f_end,f_peak,f_center,BW_3dB,BW_10dB] = freq_calc(y,y_noise,savedSignals{i,1}.fs,pt_aud,ovl);
    
    % Store all signal parameters in output structure
    signalParams.file{i,1} = savedSignals{i,1}.file; % FILENAME
    signalParams.datestr{i,1} = datestr(savedSignals{i,1}.T(1),'yyyy.mm.dd HH:MM:SS'); % DATESTRING
    signalParams.call{i,1} = savedSignals{i,1}.call; % CALL TYPE
    signalParams.comment{i,1} = savedSignals{i,1}.comment; % COMMENT
    signalParams.time{i,1} = savedSignals{i,1}.T(1); % CALL TIME (matlab datenum)
    signalParams.dur{i,1} = savedSignals{i,1}.T(end) - savedSignals{i,1}.T(1); % CALL DURATION (matlab datenum)
    signalParams.duration95{i,1} = dur/savedSignals{i,1}.fs; % 95% energy duration
    signalParams.depth{i,1} = mean(savedSignals{i,2}.p); % Depth
    signalParams.audio_peak2peak{i,1} = peak2peak(y); % Audio peak-to-peak amplitude
    signalParams.accel_mean_peak2peak{i,1} = peak2peak(A_mean); % Accel peak-to-peak amplitude (mean)
    signalParams.accel_max_peak2peak{i,1} = peak2peak(A(:,maxChannel)); % Accel peak-to-peak amplitude (max)
    signalParams.audio_rms{i,1} = rms(y); % Audio rms amplitude
    signalParams.accel_mean_rms{i,1} = rms(A_mean); % Accel rms amplitude (mean)
    signalParams.accel_max_rms{i,1} = rms(A(:,maxChannel)); % Accel rms amplitude (max)
    signalParams.audio_SNR{i,1} = y_SNR; % Audio SNR
    signalParams.accel_mean_SNR{i,1} = A_mean_SNR; % Accel SNR (mean)
    signalParams.accel_max_SNR{i,1} = A_max_SNR; % Accel SNR (max)
    signalParams.peak_frequency{i,1} = f_peak; % Peak frequency (Hz)
    signalParams.center_frequency{i,1} = f_center; % Center frequency (Hz)
    signalParams.start_frequency{i,1} = f_start; % Start frequency (Hz)
    signalParams.end_frequency{i,1} = f_end; % End frequency (Hz)
    signalParams.bandwidth3dB{i,1} = BW_3dB; % -3dB bandwidth [low, high]  <--- possibly incorrect
    signalParams.bandwidth10db{i,1} = BW_10dB;% -10dB bandwidth [low, high] <-- possibly incorrect
    
    % If saveFigs == 1, then generate figures and save as .png files
    if saveFigs
        titleStr = ['File: ',savedSignals{i,1}.file,...
            '      Time: ',datestr(savedSignals{i,1}.T(1),'yyyy.mm.dd HH:MM:SS'),...
            '      Call type: ',savedSignals{i,1}.call];
        saveStr = [saveDir, savedSignals{i,1}.call,'_',datestr(savedSignals{i,1}.T(1),'yyyy-mm-dd_HH-MM-SS'),'.png'];
        set(fpan,'Title',titleStr)
        
        waveform_MS(ax(3,1),savedSignals{i,1}.T,y);
        waveform_MS(ax(2,1),savedSignals{i,2}.T,A_mean);
        waveform_MS(ax(1,1),savedSignals{i,2}.T,A(:,maxChannel));
        
 
        freqSelection = spectrogram_MS(ax(3,2),savedSignals{i,1}.T,y,savedSignals{i,1}.fs,pt_aud,ovl,FLIM);%,[t_start t_end],[f_start f_end]);
        % WORK-AROUND TO IMPROVE START AND END FREQUENCY DETECTION
        % spectrogram_MS allows user to manually select start and end frequency
        if ~isempty(freqSelection)
            f_start = freqSelection(1); % Start frequency (Hz)
            signalParams.start_frequency{i,1} = f_start;
            f_end = freqSelection(2); % End frequency (Hz)
            signalParams.end_frequency{i,1} = f_end;
        end
        spectrogram_MS(ax(2,2),savedSignals{i,2}.T,A_mean,savedSignals{i,2}.fs,pt_acc,ovl,savedSignals{i,2}.bandpass);
        spectrogram_MS(ax(1,2),savedSignals{i,2}.T,A(:,maxChannel),savedSignals{i,2}.fs,pt_acc,ovl,savedSignals{i,2}.bandpass);
        
        powerspectrum_MS(ax(3,3),savedSignals{i,1}.T,y,savedSignals{i,1}.fs,pt_aud,ovl);
        powerspectrum_MS(ax(2,3),savedSignals{i,2}.T,A_mean,savedSignals{i,2}.fs,pt_acc,ovl,[0 savedSignals{i,1}.fs/2]);
        powerspectrum_MS(ax(1,3),savedSignals{i,2}.T,A(:,maxChannel),savedSignals{i,2}.fs,pt_acc,ovl,[0 savedSignals{i,1}.fs/2]);
        
        title(ax(3,1),'Waveform','Fontsize',8,'Fontweight','normal')
        title(ax(3,2),'Spectrogram','Fontsize',8,'Fontweight','normal')
        title(ax(3,3),'Power Spectrum','Fontsize',8,'Fontweight','normal')
        ylabel(ax(3,1),'Audio','Fontsize',8,'Fontweight','normal')
        ylabel(ax(2,1),'Accel (mean)','Fontsize',8,'Fontweight','normal')
        ylabel(ax(1,1),'Accel (max)','Fontsize',8,'Fontweight','normal')
        
        str = {['SNR_{audio} = ',num2str(y_SNR)];
            ['f_{peak} = ',num2str(f_peak)];['f_{center} = ',num2str(f_center)];
            ['f_{start} = ',num2str(f_start)];['f_{end} = ',num2str(f_end)];
            ['SNR_{accel mean} = ',num2str(A_mean_SNR)];
            ['SNR_{accel max} = ',num2str(A_max_SNR)]};
        annotate_axes(ax(2,3),str)

        set(f,'PaperPositionMode','auto')
        print(f,saveStr,'-dpng','-r300')
    end
    
end

disp('COMPLETE!')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS CALLED BY DRIVER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [SNR, dur] = snr_calc(signal,noise,varargin)
        % Modified from energy_calcs.m (Maxwell Kaplan)
        
        if isempty(varargin)
           E_trim = 0.95; % Default middle 95% of energy
        else
            E_trim = varargin{1,1}; % User input energy fraction
        end
        trim = (1-E_trim)/2;
        
        % Calculate mean-square amplitude of noise
        ms_noise = (1/length(noise))*sum(noise.^2);
        
        % Trim signal at energy bounds (default bounds 95% energy)
        E_rel = sum(signal.^2); % relative energy (sum of squared pressures)
        E_cum = cumsum(signal.^2)-ms_noise; % cumulative sum of energy (cumulative sum of squared pressures)
        E_start = find(E_cum > trim*E_rel,1);
        E_end = find(E_cum < (1-trim)*E_rel,1,'last');
        signal = signal(E_start:E_end);
        dur = length(signal); % Return the duration of the energy-bounded signal in # of samples
        
        % Calculate mean-square amplitude of signal (bounded by 95% energy)
        ms_signal = (1/length(signal))*sum(signal.^2);
        % Calculate signal-to-noise ratio in dB
        if ms_signal-ms_noise >= 0
            SNR = 10*log10((ms_signal-ms_noise)/ms_noise);
        else
            SNR = NaN;
        end
    end

    function [t_start,t_end,f_start,f_end,f_peak,f_center,BW_3dB,BW_10dB]...
            = freq_calc(signal,noise,fs,pt,ovl,varargin)
        
        if isempty(varargin)
           E_trim = 0.95; % Default middle 95% of energy
        else
            E_trim = varargin{1,1}; % User input energy fraction
        end
        trim = (1-E_trim)/2;
        
        % Calculate mean-square amplitude of noise
        ms_noise = (1/length(noise))*sum(noise.^2);
        
        % Trim signal at energy bounds (default bounds 95% energy)
        E_rel = sum(signal.^2); % relative energy (sum of squared pressures)
        E_cum = cumsum(signal.^2)-ms_noise; % cumulative sum of energy (cumulative sum of squared pressures)
        E_start = find(E_cum > trim*E_rel,1);
        E_end = find(E_cum < (1-trim)*E_rel,1,'last');
        signal = signal(E_start:E_end); % Trim signal to energy bounds
        t_start = E_start*(1/fs); % Start of energy bounds (s)
        t_end = E_end*(1/fs); % End of energy bounds (s)
        
        % Set up Fourier analysis parameters
        window = hamming(pt);
        noverlap = floor(ovl*pt);
        if length(signal) < pt % Work-around to avoid short segments causing errors
            warning('Short segment: Fourier parameters modified')
            window = hamming(length(signal));
            noverlap =  floor(ovl*length(signal));
        end
        
        % Calculate the start and end frequency from spectrogram
        [S,F,T] = spectrogram(signal,window,noverlap,pt,fs);
        [~,idx] = max(S(:,1)); % Use first column of spectrogram to determine start frequency
        f_start = mode(F(idx)); % Start frequency
        [~,idx] = max(S(:,end)); % Use last column of spectrogram to determine end frequency
        f_end = mode(F(idx)); % End frequency
        
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

    function waveform_MS(ax,t,y)
        % Function plots and formats waveform on the provided axes
        
        % Convert time vector into seconds
        t = (t-t(1))*24*60*60;
        % Plot waveform and format axis
        plot(ax,t,y)
        set(ax,'xlim',[t(1) t(end)])
    end

    function f_select = spectrogram_MS(ax,T,y,fs,pt,ovl,varargin)
        % Function plots and formats spectrogram on the provided axes
        
        % Returns manually selected start and end frequency if 95% energy
        % time bounds and automatically detected start and end frequencies
        % are provided in varargin (see below for details).
        f_select = [];
        
        % Set up Fourier analysis parameters
        window = hamming(pt);
        noverlap = floor(ovl*pt);
        if length(y) < pt % Work-around to avoid short segments causing errors
            warning('Short segment: Fourier parameters modified')
            window = hamming(length(y));
            noverlap =  floor(ovl*length(y));
        end
        
        % Calculate spectrogram
        [S,F,T] = spectrogram(y,window,noverlap,pt,fs);
        S = 10*log10(abs(S)); % Convert to decibels
        % Plot spectrogram and format axis
        imagesc(T,F,S,'Parent',ax)
        set(ax,'YDir','Normal')
        % Assign custom color map
        map = [1:-.01:0.50, 0.49:-0.03:0 ]';
        map = [map,map,map];
        colormap(ax,map)
        
       % Variable input arguments:
       % MUST be inorder: FLIM, t95, f95 -or- t95, f95 -or- FLIM
       t95 = []; f95 = [];
       switch length(varargin)
           case 1 % Only frequency axis limits provided
               if varargin{1}(1) < varargin{1}(2)
                   set(ax,'ylim',varargin{1})
               end
           case 2 % Only 95% energy time+freq bounds provided
               t95 = varargin{1};
               f95 = varargin{2};
           case 3 % Both freq axis limits and 95% energy time+freq bounds
               if varargin{1}(1) < varargin{1}(2)
                   set(ax,'ylim',varargin{1})
               end
               t95 = varargin{2};
               f95 = varargin{3}; 
       end
        
       if ~isempty(t95)
           % Plot automatically detected start and end frequencies at the
           % 95% energy time bounds (provided in varargin)
           hold(ax,'on')
           plot(ax,t95,f95,'ro');
           hold(ax,'off')
           
           % Recalculate the start frequency with manual assistance
           confirmButton = []; h = [];
           while strcmp(confirmButton,'Yes') ~= 1
               % Requires user to drag rectangle on provided axes
               delete(h); rect = getrect(ax);
               tLim = [rect(1), rect(1)+rect(3)];
               fLim = [rect(2), rect(2)+rect(4)];
               
               % Find the time indexes on the rectangle sides
               tIdxStart = find(T >= tLim(1),1,'first');
               tIdxEnd = find(T <= tLim(2),1,'last');
               % Find the freq indexes bound by the rectangle top+bottom
               fIdx = find(F >= fLim(1) & F <= fLim(2));
               
               % Calculate and store the peak frequencies on the rectangle
               % sides that fall within the frequency band specified by the
               % rectangle
               [~,fIdxStart] = max(S(fIdx,tIdxStart));
               [~,fIdxEnd] = max(S(fIdx,tIdxEnd));
               fIdxStart = fIdx(fIdxStart);
               fIdxEnd = fIdx(fIdxEnd);
               
               % Plot the rectangle and the newly selected frequencies
               hold(ax,'on')
               h(1) = rectangle('Parent',ax,'Position',rect,...
                   'EdgeColor','g','LineWidth',0.1,'LineStyle','-');
               h(2) = plot(T([tIdxStart,tIdxEnd]),F([fIdxStart,fIdxEnd]),'gx');
               hold(ax,'off')
               % Store the selected start and end frequencies
               f_select = F([fIdxStart,fIdxEnd]);
               
               % Request user to confirm selection before returning
               frStr = ['F_start =',num2str(f_select(1)),'Hz',...
                   ', F_end =',num2str(f_select(2)),' Hz'];
               confirmButton = questdlg(frStr);
           end
           
       end
    end

    function powerspectrum_MS(ax,T,y,fs,pt,ovl,varargin)
        % Function plots and formats power spectrum on the provided axes
        
        % Set up Fourier analysis parameters
        window = hamming(pt);
        noverlap = floor(ovl*pt);
        if length(y) < pt % Work-around to avoid short segments causing errors
            warning('Short segment: Fourier parameters modified')
            window = hamming(length(y));
            noverlap =  floor(ovl*length(y));
        end
        
        % Calculate Welch's periodogram
        [Pxx,Fxx]=pwelch(y,window,noverlap,pt,fs,'onesided');
        Pxx = 10*log10(abs(Pxx)); % Convert to decibels
        % Plot power spectrum
        plot(ax,Fxx,Pxx)
        
        % Suitable variable input argument provides frequency limits for
        % power spectrum
        if ~isempty(varargin)
            if varargin{1,1}(1) < varargin{1,1}(2)
                set(ax,'xlim',varargin{1})
            end
        end
    end

    function annotate_axes(ax,str)
        % Function plots string in the top right corner of provided axes
        
        ypoint = get(ax,'ylim'); ypoint = ypoint(2)-abs(0.02*ypoint(2));
        xpoint = get(ax,'xlim'); xpoint = xpoint(2)-abs(0.02*xpoint(2));
        text(xpoint,ypoint,str,'Parent',ax,...
            'HorizontalAlignment','right',...
            'VerticalAlignment','top')
    end

end