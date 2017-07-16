% decimate_for_correlation.m
% MS 2016.08.10

% Script for decimating acoustic data to the sample sample rate as
% accelerometer data

% Load-in savedSignals

for i = 1:size(savedSignals,1)
    % Get audio and accelerometer sample rates
   fs_audio = savedSignals{i,1}.fs;
   fs_accel = savedSignals{i,2}.fs;
   
   % Calculate decimation factor
   r = fs_audio/fs_accel;
   
   % Get audio and accelerometer signals (buffer with noise)
   y = [savedSignals{i,1}.noise_pre; savedSignals{i,1}.y; savedSignals{i,1}.noise_post];
   A = [savedSignals{i,2}.noise_pre; savedSignals{i,2}.y; savedSignals{i,2}.noise_post];
   y = y - mean(y);
   A(:,1) = A(:,1) - mean(A(:,1));
   A(:,2) = A(:,2) - mean(A(:,2));
   A(:,3) = A(:,3) - mean(A(:,3));
   
   % Decimate the audio signal to the sample rate of the accelerometer
   y_dec = decimate(y,r);
   A_dec = A;
   
   % Filter audio signal with high-pass filter used for accel signals
   [b,a] = butter(6,savedSignals{i,2}.bandpass(1)/(fs_accel/2),'high');
   y_dec = filter(b,a,y_dec);
   
   len = min(length(y_dec),length(A_dec));
   noiseTrim = round(0.95*length(savedSignals{i,2}.noise_pre));
   y_dec = y_dec(noiseTrim:(len-noiseTrim),1);
   A_dec = A_dec(noiseTrim:(len-noiseTrim),:);
   
   % Store the decimated signals in savedSignals
   savedSignals{i,1}.y_dec = y_dec;
   savedSignals{i,1}.A_dec = A_dec;
   savedSignals{i,1}.fs_dec = fs_accel;
end