% Load-in prh file and bm15_054a_lunge_idx_10Hz.mat
% Be careful with sample times of lunge indexes
%
% MS 2017.08.19

close all
figure;
hold all

plot_idx = 1:50:length(p);
plot(plot_idx, p(plot_idx));

plot(lunge_idx, p(lunge_idx), 'ro')

xlabel('time (samples)')
set(gca, 'ydir', 'reverse')