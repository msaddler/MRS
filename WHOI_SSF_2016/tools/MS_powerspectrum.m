    function [Pxx,Fxx] = MS_powerspectrum(ax,y,fs,pt,ovl,varargin)
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
        PxxdB = 10*log10(abs(Pxx)); % Convert to decibels
        % Plot power spectrum
        plot(ax,Fxx,PxxdB)
        
        % Variable input argument: FLIM
        if ~isempty(varargin)
            if varargin{1,1}(1) < varargin{1,1}(2)
                set(ax,'xlim',varargin{1})
            end
        end
    end