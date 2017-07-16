% Script to visualize downsweep call waveform envelopes
% MS 2016.07.24

% Load in 'savedSignals' and 'signalParams'

close all;
f = figure(1); g = figure(2);
saveDir = 'G:\Mark\data_bm\dtag_analysis\downsweep_waveform_images\';

for i=1:size(savedSignals,1)
    y = savedSignals{i,1}.y;
    y = y-mean(y);
    fs = savedSignals{i,1}.fs;
    T = savedSignals{i,1}.T;
    T = (T-T(1))*24*60*60;
    
    [b,a] = butter(1,5/(fs/2),'low');
    y_filt = filter(b,a,y.^2);
    y_filt = y_filt.^(1/2);
    
    accelSNR = signalParams.accel_max_SNR{i,1};
    audioSNR = signalParams.audio_SNR{i,1};
    
    yh = hilbert(y);
    
    figure(1)
    plot(T,y); hold on;
    plot(T,abs(yh),'r','LineWidth',.02);
    plot(T,-abs(yh),'r','LineWidth',.02);
    plot(T,y_filt,'k','LineWidth',2);
    plot(T,-y_filt,'k','LineWidth',2);
    hold off;
    
    xlim([T(1) T(end)])
    title(['Audio SNR =',num2str(audioSNR),'   Accel SNR =',num2str(accelSNR)])
    
    %     if accelSNR > 0
    %         saveStr = [saveDir,num2str(accelSNR),'acc_snr_',savedSignals{i,1}.call,'.png'];
    %     else
    %         saveStr = [saveDir,'neg',num2str(abs(accelSNR)),'acc_snr_',savedSignals{i,1}.call,'.png'];
    %     end
    %     saveas(f,saveStr)
    
    
    
    trim = (1-0.95)/2;
    E_rel = sum(y.^2); % relative energy (sum of squared pressures)
    E_cum = cumsum(y.^2); % cumulative sum of energy (cumulative sum of squared pressures)
    E_start = find(E_cum > trim*E_rel,1);
    E_end = find(E_cum < (1-trim)*E_rel,1,'last');
    y_filt = y_filt(E_start:E_end); % Trim signal to energy bounds
    t_start = E_start*(1/fs); % Start of energy bounds (s)
    t_end = E_end*(1/fs); % End of energy bounds (s)
    T = T(E_start:E_end);
    
    [maxval, idx] = max(y_filt);
    T = T - T(idx);
    
    
    if  accelSNR > 2
        color = 'r';
    else
        color = [.9 .9 .9];
    end
    
    figure(2)
    hold on;
    plot(T,y_filt,'Color',color)
    plot(T,-y_filt,'Color',color)
    hold off;
    

end