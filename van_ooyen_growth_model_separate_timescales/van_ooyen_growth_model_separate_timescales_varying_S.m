function [out] = van_ooyen_growth_model_separate_timescales_varying_S()
%{
Run the activity-dependent network growth model
Inputs:
  N: (integer) number of neurons in simulation
Returns:
  out: (structure) with fields for pos, rad, A, W, V, F
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

N = 64;

inhFract = 0.20; % Fraction of N neurons that are inhibitory
inhFlag = zeros(1, N);
[~,idx] = sort(rand(1, N));
inhFlag(idx(1:round(inhFract * N))) = 1;
inhIdx = find(inhFlag == 1); % Indexes of the inhibitory neurons
excIdx = find(inhFlag == 0); % Indexes of the excitatory neurons

dT = 0.001; % Growth model time step (slow timescale)
T = 0:dT:2.0; % Growth model time vector (slow timescale)
pos = rand(N, 2); % Neuron positions
grid_size = ceil(sqrt(N)); % GRID-POSITION NEURONS
for n = 1:N % GRID-POSITION NEURONS
   x = mod(n/grid_size, 1); % GRID-POSITION NEURONS
   y = (ceil(n/grid_size)-1) / grid_size; % GRID-POSITION NEURONS
   pos(n, :) = [x+0.025*randn(), y+0.025*randn()]; % GRID-POSITION NEURONS
end % GRID-POSITION NEURONS
rad = 0.075 * ones(length(T), N); % Initial dendritic field radii

rho = zeros(1, N); % Neuron growth rate constant
rho(excIdx) = .30;
rho(inhIdx) = .10;
epsilon = zeros(1, N); % Firing rate at which neuron growth rate = 0;
epsilon(excIdx) = 0.40; % (40)
epsilon(inhIdx) = 0.40; % (40)
beta = 0.10; % Growth rate slope parameter
H = 0.1; % (- I_saturation / E_saturation), see Van Ooyen et al. (1995)
tau_V = .001; % time constant for dV/dt

S = zeros(length(T), N, N); % Synaptic strength (time, post, pre)
% S(1, excIdx, excIdx) = 180; % S_ee (e <- e)
% S(1, excIdx, inhIdx) = 480; % S_ei (e <- i)
% S(1, inhIdx, excIdx) = 480; % S_ie (i <- e)
% S(1, inhIdx, inhIdx) = 180; % S_ii (i <- i)
for ti = 1:length(T)
    t = T(ti);
    S(ti, excIdx, excIdx) = 480 / (1 + exp((0.5 - t) / 0.15)); % S_ee (e <- e)
    S(ti, excIdx, inhIdx) = 720 / (1 + exp((0.5 - t) / 0.15)); % S_ei (e <- i)
    S(ti, inhIdx, excIdx) = 480 / (1 + exp((0.5 - t) / 0.15)); % S_ie (i <- e)
    S(ti, inhIdx, inhIdx) = 360 / (1 + exp((0.5 - t) / 0.15)); % S_ii (i <- i)
end
% for ti = 1:length(T)
%     t = T(ti);
%     S(ti, excIdx, excIdx) = 480 * (t / (t + 0.6)); % S_ee (e <- e) (480)
%     S(ti, excIdx, inhIdx) = 640 * (t / (t + 0.6)); % S_ei (e <- i) (640)
%     S(ti, inhIdx, excIdx) = 480 * (t / (t + 0.6)); % S_ie (i <- e) (480)
%     S(ti, inhIdx, inhIdx) = 360 * (t / (t + 0.6)); % S_ii (i <- i) (360)
% end

A = zeros(length(T), N, N); % Neuron area overlaps (time, post, pre)
W = zeros(length(T), N, N); % Synaptic weights
V = zeros(length(T), N); % Membrane potentials
F = zeros(length(T), N); % Neuron firing rates

% Wilson-Cowan activity model parameters
WC_params.T = [0, 100];
WC_params.dt = 0.01;
WC_params.P = 4;
WC_params.Q = 0;
WC_params.alpha_E = 2.1;
WC_params.alpha_I = 1.5;
WC_params.theta_E = 7.5;
WC_params.theta_I = 3;
WC_params.inact_E = 1;
WC_params.inact_I = 0.1;

% Vectors to store mean weights (used in WC simulation)
W_ee = zeros(length(T), 1);
W_ie = zeros(length(T), 1);
W_ei = zeros(length(T), 1);
W_ii = zeros(length(T), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Driver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ti = 1:length(T) - 1
    % Get overlap between neurons from previous step
    A(ti, :, :) = get_overlap(ti, pos, rad);
    % Scale overlaps by synaptic strength constants to get synaptic weights
    W(ti, :, :) = squeeze(S(ti, :, :)) .* squeeze(A(ti, :, :));
    
    % Calculate mean synaptic weights for Wilson-Cowan simulation
    W_ee(ti) = sum(sum(W(ti, excIdx, excIdx))) / (length(excIdx) * length(excIdx));
   	W_ie(ti) = sum(sum(W(ti, inhIdx, excIdx))) / (length(inhIdx) * length(excIdx));
    W_ei(ti) = sum(sum(W(ti, excIdx, inhIdx))) / (length(excIdx) * length(inhIdx));
    W_ii(ti) = sum(sum(W(ti, inhIdx, inhIdx))) / (length(inhIdx) * length(inhIdx));
    WC_weights.ee = W_ee(ti);
    WC_weights.ie = W_ie(ti);
    WC_weights.ei = W_ei(ti);
    WC_weights.ii = W_ii(ti);
    
    % Run Wilson-Cowan simulation to get fast-scale population activity
    [t, E, I, ~, ~, ~, ~] = WC_deterministic(WC_weights, WC_params);
    
    rand_state = rand(N, length(t));
    state = zeros(N, length(t));
    for n = 1:N
        % Calculate neuron firing rate at current time step
        if inhFlag(n) == 1
            n_state = rand_state(n, :) < I' * WC_params.dt;
        else
            n_state = rand_state(n, :) < E' * WC_params.dt;
        end
        state(n, :) = n_state;
        F(ti, n) = sum(n_state) / t(end);
        % Calculate neuron growth rate at current time step
        tmpG = growth_rate_fcn(F(ti, n), beta, epsilon(n));
        dR = dT * rho(n) * tmpG;
        rad(ti+1, n) = rad(ti, n) + dR; % Set radius for next time step
        rad(ti+1, n) = max([rad(ti+1, n), 0]);
    end
    
    %%% START: PLOT FAST TIMESCALE ACTIVITY %%%
    if mod(ti, 50) == 0
        state = state';
        figure(); 
        subplot(211)
        hold all
        plot(t, E, 'b')
        plot(t, I, 'r')
        ylim([0, 1])
        title(['time step: ', num2str(ti), ' (E = blue, I = red)'])
        
        subplot(212)
        hold all
        for n = 1:N
            plot(t, state(:,n).*(state(:,n) + n), 'k.')
        end
        drawnow;
    end
    %%% END: PLOT FAST TIMESCALE ACTIVITY %%%%
    
    for n = 1:N
        % Calculate neuron membrane potential at next time step
        W_vec = squeeze(W(ti, n, :)); % Weights of all incoming synapses
        F_vec = F(ti, :); % Firing rates of neurons synapsing on 'n'
        % Euler's method to solve for membrane potential
        tmp = -V(ti, n);
        tmp = tmp + (1 - V(ti, n)) * dot(W_vec(excIdx), F_vec(excIdx));
        tmp = tmp - (H + V(ti, n)) * dot(W_vec(inhIdx), F_vec(inhIdx));
        dV = dT * (1/tau_V) * tmp;
        V(ti+1, n) = V(ti, n) + dV;
    end
    
    % Display progress of simulation
    if mod(ti, 10) == 0
        disp([num2str(ti),'/',num2str(length(T))])
    end
end

% Build out data structure to return
out.N = N;
out.T = T;
out.pos = pos;
out.rad = rad;
out.A = A;
out.W = W;
out.V = V;
out.F = F;
out.W_ee = W_ee;
out.W_ie = W_ie;
out.W_ei = W_ei;
out.W_ii = W_ii;
out.inhFlag = inhFlag;
out.excIdx = excIdx;
out.inhIdx = inhIdx;

out.params.H = H;
out.params.tau_V = tau_V;
out.params.beta = beta;
out.params.epsilon = epsilon;
out.params.rho = rho;
out.params.S = S;

out.WC_params = WC_params;

save('out.mat', 'out')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure() % Plot radii
hold all
for n = 1:N
    if inhFlag(n)
        plot(T, rad(:, n), 'r');
    else
        plot(T, rad(:, n), 'b');
    end
end
title('Neuron Radii')

figure() % Plot overlap
hold all
A_ee = zeros(1, length(T));
A_ie = zeros(1, length(T));
A_ei = zeros(1, length(T));
A_ii = zeros(1, length(T));
for ti = 1:length(T)
    A_ee(ti) = sum(sum(A(ti, excIdx, excIdx))) / (length(excIdx) * length(excIdx));
    A_ie(ti) = sum(sum(A(ti, inhIdx, excIdx))) / (length(inhIdx) * length(excIdx));
    A_ei(ti) = sum(sum(A(ti, excIdx, inhIdx))) / (length(excIdx) * length(inhIdx));
    A_ii(ti) = sum(sum(A(ti, inhIdx, inhIdx))) / (length(inhIdx) * length(inhIdx));
end
plot(T, A_ee, 'b')
plot(T, A_ie, 'k')
plot(T, A_ei, 'r')
plot(T, A_ii, 'g')
legend('A_{ee}', 'A_{ie}', 'A_{ei}', 'A_{ii}')
title('Overlap Area')

figure() % Synaptic Weights (scaled overlap areas)
hold all
W_ee = zeros(1, length(T));
W_ie = zeros(1, length(T));
W_ei = zeros(1, length(T));
W_ii = zeros(1, length(T));
for ti = 1:length(T)
    W_ee(ti) = sum(sum(W(ti, excIdx, excIdx))) / (length(excIdx) * length(excIdx));
    W_ie(ti) = sum(sum(W(ti, inhIdx, excIdx))) / (length(inhIdx) * length(excIdx));
    W_ei(ti) = sum(sum(W(ti, excIdx, inhIdx))) / (length(excIdx) * length(inhIdx));
    W_ii(ti) = sum(sum(W(ti, inhIdx, inhIdx))) / (length(inhIdx) * length(inhIdx));
end
plot(T, W_ee, 'b')
plot(T, W_ie, 'k')
plot(T, W_ei, 'r')
plot(T, W_ii, 'g')
legend('W_{ee}', 'W_{ie}', 'W_{ei}', 'W_{ii}')
title('Synaptic Weights')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [output] = get_overlap(itrT, pos, radii)
        % Returns array of overlap areas between all neurons
        % Inputs:
        %   itrT = index of current time step
        %   pos = positions of all neurons
        %   radii = radii of all neurons at current time
        output = zeros(size(pos, 1), size(pos, 1));
        for itr1 = 1:size(pos,1)
            for itr2 = (itr1 + 1):size(pos,1)
                xdist = pos(itr2,1) - pos(itr1,1);
                ydist = pos(itr2,2) - pos(itr1,2);
                dist = (xdist^2+ydist^2)^(1/2); % Distance between centers
                r1 = radii(itrT, itr1); % Radius of neuron 1
                r2 = radii(itrT, itr2); % Radius of neuron 2
                % Determine area of overlap between neuron 1 and neuron 2
                area = r2^2*acos((dist^2 + r2^2 - r1^2)/(2*dist*r2)) + ...
                    r1^2*acos((dist^2 + r1^2 - r2^2)/(2*dist*r1)) - ...
                    (1/2)*sqrt((-1*dist+r2+r1)*(dist+r2-r1)*...
                    (dist-r2+r1)*(dist+r2+r1));
                % Store the real part of the overlap in output
                output(itr1, itr2) = real(area);
                output(itr2, itr1) = real(area);
            end
        end
    end


    function [output] = growth_rate_fcn(f, beta, epsilon)
        % Sigmoidal growth rate function (Van Ooyen et al., 1995)
        % Inputs:
        %   f : neuron firing rate
        %   beta : parameter determines steepness
        %   epsilon : parameter represents firing rate such that growth is zero
        output = 1-(2/(1+exp((epsilon-f)/beta)));
    end


end