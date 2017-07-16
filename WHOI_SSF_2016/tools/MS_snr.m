function snr = MS_snr(signal,noise,varargin)
        % Modified from energy_calcs.m (Maxwell Kaplan)
        
         % Calculate mean-square amplitude of noise
        ms_noise = (1/length(noise))*sum(noise.^2);
        
        % Variable input argument: E_trim
        % Define signal energy bounds to trim signal (i.e. middle 95%)
        if ~isempty(varargin)
           E_trim = varargin{1};
           trim = (1-E_trim)/2;
           E_rel = sum(signal.^2); % relative energy (sum of squared pressures)
           E_cum = cumsum(signal.^2)-ms_noise; % cumulative sum of energy (cumulative sum of squared pressures)
           E_start = find(E_cum > trim*E_rel,1);
           E_end = find(E_cum < (1-trim)*E_rel,1,'last');
           signal = signal(E_start:E_end);
        end 
        
        % Calculate mean-square amplitude of signal
        ms_signal = (1/length(signal))*sum(signal.^2);
        
        % Calculate signal-to-noise ratio in dB
        if ms_signal-ms_noise >= 0
            snr = 10*log10((ms_signal-ms_noise)/ms_noise);
        else
            snr = NaN;
        end
    end