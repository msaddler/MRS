% analyze_output_data_MEA.m
% MS 2017.04.07s

clear; close all

% Load-in saved data
data_filename = 'rasters_200MEA5_30min_05_08_13_aval.mat';
M = load(data_filename);
F = fieldnames(M);

state_array = M.(F{1})'; % Import raster into state_array
time_vector = [1:length(state_array)] * (1/25000); % time_vector for MEA

sum_state = sum(state_array, 2);
spike_idx = find(sum_state > 0);
spike_times = time_vector(spike_idx);
spike_counts = sum_state(spike_idx);

% Define delta_t_avg as in page 11 of Benayoun et al. (2010)
delta_t_avg = sum(spike_counts) / time_vector(end); % ??????????????????????????????????
delta_t_avg = length(spike_counts) / time_vector(end); % ???????????????????????????????
% Should delta_t_avg be defined like this instead??
delta_t_avg = time_vector(end) / sum(spike_counts);

% Directly from MarcBenayoun_Thesis_2010.pdf
% Calculate maximum number of possible avalanches given timebins
n_timebins = ceil(max(spike_times)/delta_t_avg); % Max number of avals
avalanches = zeros(1, n_timebins);
aval_times = avalanches;
spike_counter = 0;
aval_counter = 0;

% Count the number of spikes in each avalanche
ISI = diff(spike_times);
for t = 1:length(ISI)
    %spike_counter = spike_counter + 1;
    spike_counter = spike_counter + spike_counts(t+1);
    if ISI(t) >= delta_t_avg
        aval_counter = aval_counter + 1;
        avalanches(aval_counter) = spike_counter;
        aval_times(aval_counter) = spike_times(t+1);
        spike_counter = 0;
    end
end

avalanches = avalanches(1:aval_counter);
aval_times = aval_times(1:aval_counter);
time_between_aval = diff(aval_times);
log_aval = log10(avalanches);

% Decide how many size bins to create
% (e.g. 10 bins of avalanche sizes between 10 and 1000 spikes)
aval_per_bin = 100;
bins = floor(length(log_aval) / aval_per_bin);
% Use hist to get the frequency of avalanches in each size bin
[freq2, bin_center2] = hist(avalanches, bins);
freq2 = freq2/sum(freq2); % Convert counts to probability

freq = zeros(1, length(freq2));
bin_center = freq;
counter = 0;
% Remove empty bins so you don't log(0).
for t = 1:bins
    if freq2(t) ~= 0
        counter = counter + 1;
        freq(t) = log10(freq2(t));
        bin_center(t) = log10(bin_center2(t));
    end
end
freq = freq(1:counter);
bin_center = bin_center(1:counter);
[r, p_val] = corrcoef(bin_center, freq);


% Calculate slope of power law using MLE instead of naive linear fit
x_min = 10;
x_max = 10000;
counter = 0;
temp = 0;
for t = 1:length(avalanches)
    if(avalanches(t) >= x_min && avalanches(t) <= x_max)
        temp = temp + log(avalanches(t) / x_min - 0.5);
        counter = counter + 1;
    end
end
slope = -1 * (1 + counter * (1 / temp));

point = min(find((bin_center2 > x_min) .* (freq2 > 0)));
y = (freq2(point) * (bin_center2(point)) ^ (-slope)) * bin_center2 .^ (slope);

% The poiss param is the expected number of spikes in the network over
% a time of length timebin (delta_t_avg).
poiss_param = length(spike_times) * delta_t_avg * 1/max(spike_times);
% This is the theoretical size distribution for a random poisson
% process.
fit = exp(-poiss_param)*(1-exp(-poiss_param)) .^ (0:99);


% Plot results
figure

subplot(2,1,1)
hold all
plot(log10(1:100), log10(fit), 'r')
plot(bin_center, freq, 'o')
xlabel('Log10(avalanche_size)')
ylabel('Log10(probability')
legend(['theoretical size distribution for a random poisson','actual size distribution']);
% process

subplot(2,1,2) %%% Individual Neuron Activity (Raster Plot)
hold all
N = size(state_array, 2);
for n = 1:N
    s_idx = find(state_array(:, n) == 1);
    plot(time_vector(s_idx), state_array(s_idx)+n, 'k.')
end
ylabel('Neuron'); ylim([0, N])
xlabel('Time')

