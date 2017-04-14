% analyze_output_data.m
% MS 2017.04.01

clear; close all

% Load-in saved data
data_filename = 'output_data_2017-04-05_4_196neurons_randompos_delayedinhibition.mat';
M = load(data_filename);
F = fieldnames(M);
out = M.(F{1});

data = out.data;
time_vector_list = data(:, 1);
state_array_list = data(:, 2);

for i = 2:size(data, 1)
    time_vector = time_vector_list{i};
    state_array = state_array_list{i};
    
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
        spike_counter = spike_counter + 1;
        %spike_counter = spike_counter + spike_counts(t+1);%
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
%     slope = -1 * (1 + counter * (1 / temp));
%     
%     point = min(find((bin_center2 > x_min) .* (freq2 > 0)));
%     y = (freq2(point) * (bin_center2(point)) ^ (-slope)) * bin_center2 .^ (slope);
%     
%     % The poiss param is the expected number of spikes in the network over
%     % a time of length timebin (delta_t_avg).
%     poiss_param = length(spike_times) * delta_t_avg * 1/max(spike_times);
%     % This is the theoretical size distribution for a random poisson
%     % process.
%     fit = exp(-poiss_param)*(1-exp(-poiss_param)) .^ (0:99);
    
    
    N = out.N;
    excIdx = out.excIdx;
    inhIdx = out.inhIdx;
    N_E = length(excIdx);
    N_I = length(inhIdx);
    if mod(i, 5)  == 0
        figure
        subplot(3, 3, [3, 6, 9])
        hold all
        %plot(log10(1:100), log10(fit), 'r')
        plot(bin_center, freq, 'o')
        xlabel('Log10(avalanche_size)')
        ylabel('Log10(probability')
        
        subplot(3, 3, [1, 2]) %%% Population Activity
        hold all
        plot(time_vector, sum(state_array,2)/N, 'k')
        plot(time_vector, sum(state_array(:, excIdx),2)/N_E, 'b')
        plot(time_vector, sum(state_array(:, inhIdx),2)/N_I, 'r')
        legend(['All (N =',num2str(N),')'], ['E (N_E =',num2str(N_E),')'], ...
            ['I (N_I =',num2str(N_I),')'])
        ylabel('Fraction of neurons firing'); ylim([0, 1])
        
        subplot(3, 3, [4, 5, 7, 8]) %%% Individual Neuron Activity (Raster Plot)
        hold all
        for n = 1:N
            plot(time_vector, state_array(:,n).*(state_array(:,n) + n), 'k.')
        end
        ylabel('Neuron'); ylim([0, N])
        xlabel('Time')
    end
    
end

%%% RETIRED CODE %%%
%     % Temporal coarse-graining of spike times
%     coarse_spike_times = zeros(size(firing_times));
%     coarse_spike_counts = zeros(size(firing_times));
%     coarse_spike_idx = 1;
%     for ti = 1:length(firing_times)-1
%         t = firing_times(ti);
%         f = firing_counts(ti);
%         tdiff = firing_times(ti+1) - firing_times(ti);
%         if tdiff > 0.0001 % Combine spikes closer than 0.1 ms
%             coarse_spike_times(coarse_spike_idx) = t;
%             coarse_spike_counts(coarse_spike_idx) = f;
%             coarse_spike_idx = coarse_spike_idx + 1;
%         end
%     end
%     coarse_spike_times = coarse_spike_times(1:coarse_spike_idx-1);
%     coarse_spike_counts = coarse_spike_counts(1:coarse_spike_idx-1);
%     % Append the last spike time if appropriate
%     if firing_times(end) - coarse_spike_times(coarse_spike_idx-1) > 0.0001
%         coarse_spike_times = [coarse_spike_times, firing_times(end)];
%         coarse_spike_counts = [coarse_spike_counts, firing_counts(end)];
%     end
