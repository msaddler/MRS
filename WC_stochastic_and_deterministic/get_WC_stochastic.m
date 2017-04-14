function [time,state, synCount] = get_WC_stochastic(...
    W, excIdx, inhIdx, params)
% MS 2017.02.23
%
% Function to simulate Wilson-Cowan network activity stochastically at the
% level of individual neurons using Gillespie's Algorithm.
%

% General Parameters
N = length(W);
N_E = length(excIdx); % Num excitatory neurons
N_I = length(inhIdx); % Num inhibitory neurons

% Network Mean Weights (for display purposes only)
mean_W_ee = mean(mean(W(excIdx, excIdx))); % 16
mean_W_ie = mean(mean(W(inhIdx, excIdx))); % 18
mean_W_ei = mean(mean(W(excIdx, inhIdx))); % 12
mean_W_ii = mean(mean(W(inhIdx, inhIdx))); % 3
disp(['Mean weights: EE=',num2str(mean_W_ee),', IE=', num2str(mean_W_ie),...
    ', EI=',num2str(mean_W_ei),', II=', num2str(mean_W_ii)])

% Model Parameters
T = params.T(2); % Epoch length (s)

h = zeros(1,N); % Input specified for each neuron
h(excIdx) = params.P;
h(inhIdx) = params.Q;

inact = zeros(1,N); % Decay probability of the active state of the neuron
inact(excIdx) = params.inact_E;
inact(inhIdx) = params.inact_I;

alpha = zeros(1,N); % Firing rate function parameter: steepness
theta = zeros(1,N); % Firing rate function parameter: threshold
alpha(excIdx) = params.alpha_E;
alpha(inhIdx) = params.alpha_I;
theta(excIdx) = params.theta_E;
theta(inhIdx) = params.theta_I;

count = 1; % Counter for stochastic arrays
cum_t = 0; % Initial time
time(count) = cum_t; % Vector for stochastic timebase
state(count, 1:N) = rand(1,N) < 0.00; % Array for neuron states (row = time, col = cell)

while cum_t < T
    % Add next row to state array
    if count > 1
        state(count, :) = state(count - 1, :);
    end
    
    % Compute the neuronal transition rates
    tr = zeros(1, N);
    for n = 1:N
        if state(count, n) == 1 % if neuron is active
            tr(1, n) = inact(n);
        else % if neuron is inactive
            s = (1 / N_E)*(dot(W(n, excIdx), state(count, excIdx)))... % Add excitatory inputs
                - (1 / N_I)*(dot(W(n, inhIdx), state(count, inhIdx)))... % Subtract inhibitory inputs
                + h(n); % Add external inputs
            %f = firing_rate_response(s);
            %f = sigmoid_firing_rate(s, alpha(n), theta(n));
            f = sigmoid_firing_rate(s, alpha(n), theta(n));
            tr(1, n) = f;
        end
    end
    
    % Compute network transition rate and cumulative transition rates
    network_tr = sum(tr); % calculate network transition rate
    tr = tr ./ network_tr; % normalize neuronal transition rates
    cumulative_tr = cumsum(tr); % cumulative function of all probabilities
                                % associated with the cells in N
    
    % Pick random numbers from a uniform distribution
    r = rand(1, 2);
    
    % Pick time increment (dt) from an exponential distribution of rate ntr
    dt = (1 / network_tr) * log(1 / r(1)); % Gillespie (1977)
    
    % Now pick a neuron (mu) to update using the cumulative distribution
    pick = 0;
    pick_index = 0;
    while pick == 0
        pick_index = pick_index + 1;
        if cumulative_tr(pick_index) >= r(2) % Gillespie (1977)
            mu = pick_index;
            pick = 1;
        end
    end
    
    % Update the state of neuron mu and move to next time step
    state(count, mu) = mod(state(count, mu) + 1, 2);
    time(count) = cum_t;
    cum_t = cum_t + dt;
    count = count + 1;
    
end
disp('Gillespie algorithm complete')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

preExcIdx = excIdx; %find(sum(state(:, excIdx)) > 0);
exc_synapses = W(:, preExcIdx);
exc_synapses = exc_synapses(exc_synapses > 0);
exc_count = size(exc_synapses, 1);

preInhIdx = inhIdx; %find(sum(state(:, inhIdx)) > 0);
inh_synapses = W(:, preInhIdx);
inh_synapses = inh_synapses(inh_synapses > 0);
inh_count = size(inh_synapses, 1);

synCount.E = exc_count;
synCount.I = inh_count;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function f = firing_rate_response(s)
        % f is the response function, giving the firing rate as a function
        % of input, s (Benayoun et al., 2010)
        %
        % Note that we use a threshold of 4.5 for the input (to align
        % better with Meijer et al. (2015)
        if s > 4.5
            f = tanh(s-4.5);
        else
            f = 0;
        end
    end

    function [output] = sigmoid_firing_rate(X,alpha,theta)
        % Sigmoid firing rate function (Wilson and Cowan, 1972)
        %   Parameters set location and value of maximum slope:
        %   max[S'(X)] = S'(theta) = alpha/4
        
        output = (1/(1+exp(-alpha*(X-theta)))) - (1/(1+exp(alpha*theta)));
        %output = 1 / (1 + exp((theta - X) / alpha)); %(Van Ooyen 1994)
    end

    function [output] = gaussian_firing_rate(Jxk,Xsd,Xtheta)
        % Gaussian firing rate function (Meijer et al., 2015)
        output = exp(-1*((Jxk-Xtheta)/Xsd)^2)-exp(-(-Xtheta/Xsd)^2);
    end

end