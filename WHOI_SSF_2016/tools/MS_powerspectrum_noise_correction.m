function [Pxx,Fxx,Pyy,Fyy] = MS_powerspectrum_noise_correction(ax,x,y,fs,pt,ovl,varargin)
% Function plots and formats power spectrum on the provided axes
% x := signal, y := noise

% Set up Fourier analysis parameters
window = hamming(pt);
noverlap = floor(ovl*pt);
if length(x) < pt % Work-around to avoid short segments causing errors
    warning('Short segment: Fourier parameters modified')
    window = hamming(length(x));
    noverlap =  floor(ovl*length(x));
end

% Calculate Welch's periodogram of signal
[Pxx,Fxx]=pwelch(x,window,noverlap,pt,fs,'onesided');

%%%%%%%%%%%%%%%%%%%%
% NOISE CORRECTION %
%%%%%%%%%%%%%%%%%%%%

if length(y) < pt % Work-around to avoid short segments causing errors
    warning('Short segment: Fourier parameters modified')
    window = hamming(length(y));
    noverlap =  floor(ovl*length(y));
end

% Calculate Welch's periodogram of noise
[Pyy,Fyy]=pwelch(y,window,noverlap,pt,fs,'onesided');

%%%%%%%%%%%%%%%%%%%%
% NOISE CORRECTION %
%%%%%%%%%%%%%%%%%%%%

% Return power spectrum of signal - power spectrum of noise
%Pxx = Pxx - Pyy;

% Plot power spectrum
PxxdB = 10*log10(abs(Pxx)); % Convert to decibels
PyydB = 10*log10(abs(Pyy)); % Convert to decibels
plot(ax,Fxx,PxxdB); hold(ax,'on') % Plot signal spectrum
plot(ax,Fyy,PyydB,'r-.'); hold(ax,'off') % Plot noise spectrum
%legend(ax,'Signal','Noise','Location','Southeast')

% Variable input argument: FLIM
if ~isempty(varargin)
    if varargin{1,1}(1) < varargin{1,1}(2)
        set(ax,'xlim',varargin{1})
    end
end

end