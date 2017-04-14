function out = growth_model_WC_stochastic_and_deterministic()
close all

% Network growth model time vector
dT = 0.01; % Time step
T = (0:dT:1)'; % Time vector

% Initialize output data
out.data = cell(length(T), 2); % Each row: time_gillespie, state

% Assign number and distribution of excitatory/inhibitory neurons
N = 196;
inhFract = 0.20; % Fraction of N neurons that are inhibitory
inhFlag = zeros(1, N);
[~,idx] = sort(rand(1, N));
inhFlag(idx(1:round(inhFract * N))) = 1;
excIdx = find(inhFlag == 0); % Indexes of the excitatory neurons
inhIdx = find(inhFlag == 1); % Indexes of the inhibitory neurons

% Initialize arrays and vectors
pos = rand(N,2); % Neuron positions
% grid_size = ceil(sqrt(N));
% for n = 1:N % GRID-POSITION NEURONS
%    x = mod(n/grid_size, 1);
%    y = (ceil(n/grid_size)-1) / grid_size;
%    pos(n, :) = [x+0.025*randn(), y+0.025*randn()];
% end
rad = 0.01 * ones(length(T),N); % Dendritic field radii
A = zeros(length(T), N, N); % Neuron area overlaps (time,pre,post)
W = zeros(length(T), N, N); % Synaptic weights
F = zeros(length(T), N); % Membrane potentials

S = zeros(N, N); % Synaptic strength constants (row = post, col = pre)
S(excIdx, excIdx) = 320000 * 2.5 / 100; % S_ee : E <- E
S(inhIdx, excIdx) = 360000 * 2.5 / 100; % S_ie : I <- E
S(excIdx, inhIdx) = 240000 * 2.5 / 100; % S_ei : E <- I
S(inhIdx, inhIdx) = 60000 * 2.5 / 100; % S_ii : I <- I
S = S / (N); % Changed from N^2

% Network growth model parameters
rho = zeros(1, N);
rho(excIdx) = 1.00/0.5; % Define growth rate const for excitatory neurons
rho(inhIdx) = 1.00/0.5; % Define growth rate const for inhibitory neurons
beta = 0.20; % Growth rate slope parameter
epsilon = zeros(1, N);
epsilon(excIdx) = 0.35; % Firing rate at which growth rate = 0
epsilon(inhIdx) = 0.35; % (Quite arbitrarily chosen!!!)

% Wilson-Cowan network activity parameters
WC_params.T = [0 300];
WC_params.dt = 0.01;
WC_params.P = 4.00;
WC_params.Q = 0.00;
WC_params.alpha_E = 2*1.5828; % Meijer et al. (2015), sigmoid firing rate
WC_params.alpha_I = 2*2.2201;
WC_params.theta_E = 5.2516 * 1.1;
WC_params.theta_I = 3.7512 * 0.9;
WC_params.inact_E = .25;
WC_params.inact_I = .05;

% Realtime trace figure
fig = figure('Position', [50 50 500 900]);
ax_handles(1) = subplot(3,1,1,'Parent', fig);
ax_handles(2) = subplot(3,1,2,'Parent', fig);
ax_handles(3) = subplot(3,1,3,'Parent', fig);
hold(ax_handles(1), 'on'); ylabel(ax_handles(1), 'synapses')
hold(ax_handles(2), 'on'); ylabel(ax_handles(2), 'radii')
hold(ax_handles(3), 'on'); ylabel(ax_handles(3), 'weights')
set(ax_handles(1), 'xlim', [T(1) T(end)])
set(ax_handles(2), 'xlim', [T(1) T(end)])
set(ax_handles(3), 'xlim', [T(1) T(end)])

for i = 2:length(T)
    % Get overlap between neurons from previous step
    A(i,:,:) = get_overlap(pos, rad(i-1, :));
    % Scale overlaps by synaptic strength constants to get synaptic weights
    Wi = S .* squeeze(A(i, :, :));
    W(i,:,:) = Wi;
    
    % Get individual neuron activity using Gillespie's Algorithm
    [time_gillespie, state, synCount] = get_WC_stochastic(...
        Wi, excIdx, inhIdx, WC_params);
    out.data{i, 1} = time_gillespie;
    out.data{i, 2} = state;
    
    % Get deterministic population activity using Wilson-Cowan equations
    [time_deterministic, E, I] = get_WC_deterministic(...
        Wi, excIdx, inhIdx, WC_params);
    
    % Calculate neuron firing rate and update growth model
    %F(i,:) = sum(state)/length(time_gillespie); % (used 2017-04-05)
    F(i,:) = sum(state)/time_gillespie(end); % (used 2017-04-10)
    
    disp(['Mean Network Activity = ', num2str(mean(F(i,:)))])
    for n = 1:N
        G = get_growth_rate(F(i, n), beta, epsilon(n));
        dR = dT * rho(n) * G;
        rad(i,n) = max(rad(i-1, n) + dR, 0);
    end
    
    % Generate plots!
    if mod(T(i), 0.05) == 0
        growth_model_WC_plot_generator(time_deterministic, E, I, ...
            time_gillespie, state, excIdx, inhIdx, T(i))
        weights.ee = mean(mean(Wi(excIdx, excIdx))); % 16
        weights.ie = mean(mean(Wi(inhIdx, excIdx))); % 18
        weights.ei = mean(mean(Wi(excIdx, inhIdx))); % 12
        weights.ii = mean(mean(Wi(inhIdx, inhIdx))); % 3
        f = get_WC_phase_plane(weights, WC_params);
        E_stoch = sum(state(:, excIdx), 2) / length(excIdx);
        I_stoch = sum(state(:, inhIdx), 2) / length(inhIdx);
        plot(E_stoch, I_stoch)
    end
    % Realtime trace plotting
    growth_model_WC_plot_realtime(ax_handles, T(i), synCount, rad(i,:), ...
        Wi, excIdx, inhIdx)
    
end

out.N = N;
out.excIdx = excIdx;
out.inhIdx = inhIdx;
out.pos = pos;
out.rad = rad;
out.T = T;
out.dT = dT;
out.S = S;
out.A = A;
out.W = W;
out.F = F;
out.rho = rho;
out.epsilon = epsilon;
out.beta = beta;
out.WC_params = WC_params;

% SAVE OUTPUT VARIABLE
save('TEMP.mat', 'out')

figure; hold all
for i = 1:out.N
    if inhFlag(i)
        plot(out.T, out.rad(:, i), 'r')
    else
        plot(out.T, out.rad(:, i), 'k')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function output = get_overlap(positions,radii)
        % Function calculate the fractional overlaps between all pairs of
        % neurons. The fractional overlap between a pre-synaptic (pre) and
        % post-synaptic (post) neuron is calculated as:
        %   output(post, pre) = circ_overlap(post, pre) / circ_area(pre)
        %   output(pre, post) = circ_overlap(post, pre) / circ_area(post)
        %
        % INPUTS:
        %   positions: positions of all neurons (N-by-2 array)
        %   radii: radii of all neurons (vector of length N)
        % RETURNS:
        %   output: matrix of fractional overlaps (row = post, col = pre)
        
        output = zeros(size(positions, 1), size(positions, 1));
        for row = 1:size(positions, 1)
            for col = (row + 1):size(positions, 1)
                xdist = positions(col,1) - positions(row,1);
                ydist = positions(col,2) - positions(row,2);
                dist = (xdist^2+ydist^2)^(1/2); % Distance between centers
                ra = radii(row); % Radius of neuron "a" (row = pre)
                rb = radii(col); % Radius of neuron "b" (col = post)
                % Determine area of overlap between pre and post circle
                if ra <= 0 || rb <= 0
                    area = 0;
                else
                    area = rb^2*acos((dist^2 + rb^2 - ra^2)/(2*dist*rb)) + ...
                        ra^2*acos((dist^2 + ra^2 - rb^2)/(2*dist*ra)) - ...
                        (1/2)*sqrt((-1*dist+rb+ra)*(dist+rb-ra) * ...
                        (dist-rb+ra)*(dist+rb+ra));
                end
                % Store the real part of the overlap in output
                output(row, col) = real(area) / (pi * rb^2);
                output(col, row) = real(area) / (pi * ra^2);
            end
        end
    end


    function output = get_growth_rate(f, beta, epsilon)
        % Sigmoidal growth rate function (Van Ooyen et al., 1995)
        % MS 2016.08.23
        
        % f : neuron firing rate
        % beta : parameter determines steepness
        % epsilon : parameter represents firing rate s.t. growth is zero
        
        output = 1 - (2 / (1 +exp((epsilon - f) / beta)));
    end


end