function [S,F,T] = MS_spectrogram(ax,y,fs,pt,ovl,varargin)
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
colormap('Bone') % Default colormap

% Variable input arguments:
% FLIM (Frequency limits for spectrogram), MAP (custom color map)
switch length(varargin)
    case 1 % FLIM provided
        if varargin{1}(1) < varargin{1}(2)
            set(ax,'ylim',varargin{1})
        end
    case 2 % FLIM and MAP provided
        if varargin{1}(1) < varargin{1}(2)
            set(ax,'ylim',varargin{1})
        end
        map = varargin{2};
        colormap(ax,map)
end

end