function [signalParams,signalContours] = dtaganalyze_contours(savedSignals)
% MS 2016.07.27
% Improved version of dtaganalyze.m for calculating peak frequency contours
% along tonal calls (blue whale downsweeps) in both audio and accelerometer
% data.
%
% Script for calculating signal parameters from combined audio and
% accerlerometer signals saved using dtagsignal.m

% Initialize output variables
signalParams = struct();
signalContours = struct();

% Fourier Analysis Parameters
pt_aud = 256; % FFT size (audio data)
pt_acc = 128; % FFT size (accelerometer data)
ovl = 0.99; % FFT overlap
FLIM = [0 500]; % Limits for audio spectrogram

saveDir = 'G:\Mark\data_bm\dtag_analysis\dtaganalyze_images\';
tmpSaveStr = [cd,'\TEMP.mat'];

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
    [t95_start,t95_end,~,~,f_peak,f_center,BW_3dB,BW_10dB] = freq_calc(y,y_noise,savedSignals{i,1}.fs,pt_aud,ovl);
    % Calculate frequency parameters for accel signal
    [~,~,~,~,f_peak_acc,f_center_acc,~,~] = freq_calc(A(:,maxChannel), A_noise(:,maxChannel),savedSignals{i,2}.fs,pt_acc,ovl);
    
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
    signalParams.audio_f_peak{i,1} = f_peak; % Audio Peak frequency (Hz)
    signalParams.audio_f_center{i,1} = f_center; % Audio Center frequency (Hz)
    signalParams.audio_f_start{i,1} = NaN; % Audio Start frequency (Hz)             <-- NaN for now
    signalParams.audio_f_end{i,1} = NaN; % Audio End frequency (Hz)                 <-- NaN for now
    signalParams.accel_f_peak{i,1} = f_peak_acc; % Accel Peak Frequency (Hz)
    signalParams.accel_f_center{i,1} = f_center_acc; % Accel Center Frequency (Hz)
    signalParams.accel_f_start{i,1} = NaN; % Accel Start Frequency (Hz)             <-- NaN for now
    signalParams.accel_f_end{i,1} = NaN; % Accel End Frequency (Hz)                 <-- NaN for now
    signalParams.bandwidth3dB{i,1} = BW_3dB; % -3dB bandwidth [low, high]           <--- possibly incorrect
    signalParams.bandwidth10db{i,1} = BW_10dB;% -10dB bandwidth [low, high]         <-- possibly incorrect
    
    % Generate plots to save as .png file
    titleStr = ['File: ',savedSignals{i,1}.file,...
        '      Time: ',datestr(savedSignals{i,1}.T(1),'yyyy.mm.dd HH:MM:SS'),...
        '      Call type: ',savedSignals{i,1}.call];
    saveStr = [saveDir,'Contours_',savedSignals{i,1}.call,'_',datestr(savedSignals{i,1}.T(1),'yyyy-mm-dd_HH-MM-SS'),'.png'];
    set(fpan,'Title',titleStr)
    
    % Plot all signal waveforms
    waveform_MS(ax(3,1),savedSignals{i,1}.T,y);
    waveform_MS(ax(2,1),savedSignals{i,2}.T,A_mean);
    waveform_MS(ax(1,1),savedSignals{i,2}.T,A(:,maxChannel));
    % Plot Audio Signal Spectrogram
    TMARKS = [(t95_start+t95_end)/2, (t95_start+t95_end)/2];
    FMARKS = [f_peak, f_center];
    [S_aud,F_aud,T_aud] = spectrogram_MS(ax(3,2),y,savedSignals{i,1}.fs,pt_aud,ovl,...
        FLIM,TMARKS,FMARKS);
    % Plot Mean Accel Signal Spectrogram
    spectrogram_MS(ax(2,2),A_mean,savedSignals{i,2}.fs,pt_acc,ovl,...
        savedSignals{i,2}.bandpass);
    % Plot Max Accel Signal Spectrogram
    TMARKS = [(t95_start+t95_end)/2, (t95_start+t95_end)/2];
    FMARKS = [f_peak_acc, f_center_acc];
    [S_acc,F_acc,T_acc] = spectrogram_MS(ax(1,2),A(:,maxChannel),savedSignals{i,2}.fs,pt_acc,ovl,...
        savedSignals{i,2}.bandpass,TMARKS,FMARKS);
    % Plot all signal power spectra
    powerspectrum_MS(ax(3,3),y,savedSignals{i,1}.fs,pt_aud,ovl);
    powerspectrum_MS(ax(2,3),A_mean,savedSignals{i,2}.fs,pt_acc,ovl,[0 savedSignals{i,1}.fs/2]);
    powerspectrum_MS(ax(1,3),A(:,maxChannel),savedSignals{i,2}.fs,pt_acc,ovl,[0 savedSignals{i,1}.fs/2]);
    
    % Format plots
    title(ax(3,1),'Waveform','Fontsize',8,'Fontweight','normal')
    title(ax(3,2),'Spectrogram','Fontsize',8,'Fontweight','normal')
    title(ax(3,3),'Power Spectrum','Fontsize',8,'Fontweight','normal')
    ylabel(ax(3,1),'Audio','Fontsize',8,'Fontweight','normal')
    ylabel(ax(2,1),'Accel (mean)','Fontsize',8,'Fontweight','normal')
    ylabel(ax(1,1),'Accel (max)','Fontsize',8,'Fontweight','normal')
    
    % Request user input in finding peak frequency contours for audio and
    % max accel signals
    [t_audio_contour,f_audio_contour] = get_contour(ax(3,2),S_aud,F_aud,T_aud);
    [t_accel_contour,f_accel_contour] = get_contour(ax(1,2),S_acc,F_acc,T_acc);
    
    f_start_aud = [];
    f_end_aud = [];
    f_start_acc = [];
    f_end_acc = [];
    
    % Store audio contour in output data structure and update the audio
    % start and end frequencies
    if ~isempty(f_audio_contour)
        f_start_aud = f_audio_contour(1);
        f_end_aud = f_audio_contour(end);
        signalParams.audio_f_start{i,1} = f_start_aud; % Audio Start frequency (Hz)
        signalParams.audio_f_end{i,1} = f_end_aud; % Audio End frequency (Hz)
        signalContours.t_audio{i,1} = t_audio_contour;
        signalContours.f_audio{i,1} = f_audio_contour;
    end
    % Store accel contour in output data structure and update the accel
    % start and end frequencies
    if ~isempty(f_accel_contour)
        f_start_acc = f_accel_contour(1);
        f_end_acc = f_accel_contour(end);
        signalParams.accel_f_start{i,1} = f_start_acc; % Accel Start Frequency (Hz)
        signalParams.accel_f_end{i,1} = f_end_acc; % Accel End Frequency (Hz)
        signalContours.t_accel{i,1} = t_accel_contour;
        signalContours.f_accel{i,1} = f_accel_contour;
    end
    
    % Annotate axes with call parameters
    str = {
        ['f_{peak,audio} = ',num2str(f_peak),'   ','f_{center,audio} = ',num2str(f_center)];
        ['f_{peak,accel} = ',num2str(f_peak_acc),'   ','f_{center,accel} = ',num2str(f_center_acc)];
        ['f_{start,audio} = ',num2str(f_start_aud),'   ','f_{end,audio} = ',num2str(f_end_aud)];
        ['f_{start,accel} = ',num2str(f_start_acc),'   ','f_{end,accel} = ',num2str(f_end_acc)];
        ['SNR_{audio} = ',num2str(y_SNR)];
        ['SNR_{accel mean} = ',num2str(A_mean_SNR)];
        ['SNR_{accel max} = ',num2str(A_max_SNR)]
        };
    annotate_axes(ax(2,3),str)
    
    % Save plots as .png image
    set(f,'PaperPositionMode','auto')
    print(f,saveStr,'-dpng','-r300')
    
    % Save output structures at each step to minimize lost time
    save(tmpSaveStr,'signalParams','signalContours');
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
        [S,F,~] = spectrogram(signal,window,noverlap,pt,fs);
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

    function [S,F,T] = spectrogram_MS(ax,y,fs,pt,ovl,varargin)
        % Function plots and formats spectrogram on the provided axes
        
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
        % Note that t95 and f95 can be any pair of two points to plot on
        % top of the spectrogram (i.e. center and peak frequency)
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
            % Plot the provided pair of frequencies of interest
            hold(ax,'on')
            plot(ax,t95,f95,'rx');
            hold(ax,'off')
        end
    end

    function powerspectrum_MS(ax,y,fs,pt,ovl,varargin)
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

    function [t_contour,f_contour] = get_contour(ax,S,F,T)
        confirmButton = []; h = [];
        while strcmp(confirmButton,'Yes') ~= 1
            % Requires user to drag rectangle on provided axes
            
            if strcmp(confirmButton,'Cancel') == 1
                t_contour = [];
                f_contour = [];
                return
            end
            
            delete(h); rect = getrect(ax);
            tLim = [rect(1), rect(1)+rect(3)];
            fLim = [rect(2), rect(2)+rect(4)];
            
            % Find the time indexes bound by rectangle sides
            tIdx = find(T >= tLim(1) & T <= tLim(2));
            % Find the freq indexes bound by the rectangle top+bottom
            fIdx = find(F >= fLim(1) & F <= fLim(2));
            
            % Calculate and store the peak frequencies on the rectangle
            % sides that fall within the frequency band specified by the
            % rectangle
            [~,fPeakIdx] = max(S(fIdx,tIdx));
            fPeakIdx = fIdx(fPeakIdx);
            
            t_contour = T(tIdx);
            f_contour = F(fPeakIdx);
            
            % Plot the rectangle and the newly selected frequencies
            try
                hold(ax,'on')
                h(1) = rectangle('Parent',ax,'Position',rect,...
                    'EdgeColor','g','LineWidth',0.1);
                h(2) = plot(t_contour,f_contour,'go');
                hold(ax,'off')
            catch
                warning('Contour plotting error')
            end
            
            % Request user to confirm selection before returning
            frStr = ['F_start =',num2str(f_contour(1)),'Hz',...
                ', F_end =',num2str(f_contour(end)),' Hz'];
            confirmButton = questdlg(frStr);
        end
    end

end