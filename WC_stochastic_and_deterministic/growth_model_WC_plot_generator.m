function growth_model_WC_plot_generator(t_det, E_det, I_det, ...
                           time, state, excIdx, inhIdx, varargin)

N = size(state, 2);
N_E = length(excIdx);
N_I = length(inhIdx);

figure;

subplot(4,1,1) %%% Deterministic Model
hold all
plot(t_det, E_det, 'b')
plot(t_det, I_det, 'r')
legend('E', 'I')
ylabel('Population Activity'); ylim([0, 1])
title(['Gillespie Network Activity: N = ', num2str(N), ...
       ', T = ', num2str(varargin{1})])

subplot(4,1,2:3) %%% Population Activity
hold all
plot(time, sum(state,2)/N, 'k')
plot(time, sum(state(:, excIdx),2)/N_E, 'b')
plot(time, sum(state(:, inhIdx),2)/N_I, 'r')
legend(['All (N =',num2str(N),')'], ['E (N_E =',num2str(N_E),')'], ...
    ['I (N_I =',num2str(N_I),')'])
ylabel('Fraction of neurons firing'); ylim([0, 1])

subplot(4,1,4) %%% Individual Neuron Activity (Raster Plot)
hold all
for n = 1:N
    plot(time, state(:,n).*(state(:,n) + n), 'k.')
end
ylabel('Neuron'); ylim([0, N])
xlabel('Time')
drawnow;

end