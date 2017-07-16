clear; close all

load('/Users/mark/Documents/WHOI_SSF_2016/Mark/data_bm/dtag_signals/bm15_054a_DS_calls_v3_noise_1s.mat')

aud_start_times = [];
acc_start_times = [];
for i = 1:5
    aud_vec = AUD.T{i};
    acc_vec = ACC.T{i};
    aud_start_times = [aud_start_times, aud_vec(1)];
    acc_start_times = [acc_start_times, acc_vec(1)];
end

for i = 1:size(savedSignals, 1)
    
    audioSignal = savedSignals{i, 1};
    accelSignal = savedSignals{i, 2};
    
    audVecIdx = find(aud_start_times < audioSignal.T(1), 1, 'last');
    accVecIdx = find(acc_start_times < accelSignal.T(1), 1, 'last');
    
    aud_vec = AUD.T{audVecIdx};
    acc_vec = ACC.T{accVecIdx};
    aud_y = AUD.y{audVecIdx};
    acc_y = ACC.y{accVecIdx};
    
    aud_n = length(audioSignal.y); % # of samples to use for audio noise
    acc_n = length(accelSignal.y); % # of samples to use for accel noise
    
    audPreIdx = find(aud_vec < audioSignal.T(1), aud_n, 'last');
    accPreIdx = find(acc_vec < accelSignal.T(1), acc_n, 'last');
    audPostIdx = find(aud_vec > audioSignal.T(length(audioSignal.T)), aud_n);
    accPostIdx = find(acc_vec > accelSignal.T(length(accelSignal.T)), acc_n);
    
    % Modify noise_pre and noise_post with the new noise indexes
    audioSignal.noise_pre = aud_y(audPreIdx, 1);
    audioSignal.noise_post = aud_y(audPostIdx, 1);
    accelSignal.noise_pre = acc_y(accPreIdx, :);
    accelSignal.noise_post = acc_y(accPostIdx, :);
    
    % Update savedSignals with newly modified audioSignal and accelSignal
    savedSignals{i, 1} = audioSignal;
    savedSignals{i, 2} = accelSignal;
    
    disp(i)
end