signal_p2p = [];

mean_p2p_noise_pre = [];
mean_p2p_noise_post = [];
mean_p2p_noise = [];

for i = 1:size(savedSignals, 1)

    aud_sig = savedSignals{i, 1};
    acc_sig = savedSignals{i, 2};
    
    aud = aud_sig.y;
    acc = acc_sig.y;
    acc_p2p = peak2peak(acc);
    signal_p2p = [signal_p2p; mean(acc_p2p)];
    
    acc_noise_pre = acc_sig.noise_pre;
    acc_noise_post = acc_sig.noise_post;
    
    p2p_pre = peak2peak(acc_noise_pre);
    p2p_post = peak2peak(acc_noise_post);
    p2p_noise = mean([p2p_pre; p2p_post]);
    
    mean_p2p_noise_pre = [mean_p2p_noise_pre; (p2p_pre)];
    mean_p2p_noise_post = [mean_p2p_noise_post; (p2p_post)];
    mean_p2p_noise = [mean_p2p_noise; mean(p2p_noise)];
end

disp(mean(mean_p2p_noise_pre))
disp(std(mean_p2p_noise_pre))

disp(mean(mean_p2p_noise_post))
disp(std(mean_p2p_noise_post))

figure; hold all
plot(mean_p2p_noise_pre(:, 1), 'r.')
plot(mean_p2p_noise_pre(:, 2), 'g.')
plot(mean_p2p_noise_pre(:, 3), 'b.')
plot(mean_p2p_noise_post(:, 1), 'r.')
plot(mean_p2p_noise_post(:, 2), 'g.')
plot(mean_p2p_noise_post(:, 3), 'b.')

disp(mean(mean_p2p_noise))
disp(std(mean_p2p_noise))
disp([min(mean_p2p_noise), max(mean_p2p_noise)])

figure; hold all
plot(mean_p2p_noise, 'o')
plot(signal_p2p, 'x')

disp(mean(signal_p2p))
disp(std(signal_p2p))