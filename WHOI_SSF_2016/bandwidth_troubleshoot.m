% Checking -3dB and -10dB bandwidths (2017.07.22)

close all
num = 61;

signal = savedSignals{num,1}.y;
signal = signal-mean(signal);
noise = savedSignals{num,1}.noise_pre;
noise = noise-mean(noise);
fs = savedSignals{num,1}.fs;

pt = 256;
ovl = 0.99;

E_trim = 0.95;
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

% Calculate the start and end frequency from spectrogram
[S,F,T] = spectrogram(signal,hamming(pt),floor(ovl*pt),pt,fs);
[~,idx] = max(S(:,1)); % Use first column of spectrogram to determine start frequency
f_start = mode(F(idx)); % Start frequency
[~,idx] = max(S(:,end)); % Use last column of spectrogram to determine end frequency
f_end = mode(F(idx)); % End frequency

% Calculate peak and center frequencies from periodogram
[Pxx,Fxx]=pwelch(signal,hamming(pt),floor(ovl*pt),pt,fs,'onesided');
[~,peakIdx] = max(Pxx); % Find where power spectrum has maximum
f_peak = Fxx(peakIdx); % Peak frequency
f_center = sum(Pxx.*Fxx)/sum(Pxx); % Center frequency

% Calculate -3dB and -10db bandwidths
PxxdB = 10*log10(Pxx); % Convert power spectrum to dB
ampPeak = PxxdB(peakIdx); % Peak amplitude in dB

BW = -3;
idx3 = find(PxxdB-ampPeak > BW);
db3_BW = [min(Fxx(idx)), max(Fxx(idx))];

BW = -10;
idx10 = find(PxxdB-ampPeak > BW);
db10_BW = [min(Fxx(idx)), max(Fxx(idx))];

% From Max Kaplan's whistle_analysis_MBK_DEC2014.m
pAmpl = max(PxxdB); % What is the peak amplitude
[Rpeak,~] = find(PxxdB==max(PxxdB)); % What is the row index for the peak amplitude?

% Let's do the -3dB BW first
BW = 3;
low3 = find(abs(PxxdB(1:Rpeak)-(pAmpl-BW))==min(abs(PxxdB(1:Rpeak)-(pAmpl-BW))));
high3 = Rpeak+find(abs(PxxdB(Rpeak:end)-(pAmpl-BW))==min(abs(PxxdB(Rpeak:end)-(pAmpl-BW))));
BW3DB(1,1) = F(low3);
BW3DB(1,2) = F(high3);

% And now for the -10dB BW:
BW = 10;
low10 = find(abs(PxxdB(1:Rpeak)-(pAmpl-BW))==min(abs(PxxdB(1:Rpeak)-(pAmpl-BW))));
high10 = Rpeak+find(abs(PxxdB(Rpeak:end)-(pAmpl-BW))==min(abs(PxxdB(Rpeak:end)-(pAmpl-BW))));
BW10DB(1,1) = F(low10);
BW10DB(1,2) = F(high10);

figure; hold all
plot(Fxx,PxxdB,'.')
plot((pAmpl-3)*ones(size(Fxx)),'r')
plot((pAmpl-10)*ones(size(Fxx)),'b')

plot(Fxx(idx3([1 end])),PxxdB(idx3([1 end])),'rx')
plot(Fxx(idx10([1 end])),PxxdB(idx10([1 end])),'bx')

plot(BW3DB,[PxxdB(low3) PxxdB(high3)],'ro')
plot(BW10DB,[PxxdB(low10) PxxdB(high10)],'bo')