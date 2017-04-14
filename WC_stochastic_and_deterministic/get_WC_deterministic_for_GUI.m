function [t, E, I, J_E, J_I, F_E, F_I] = get_WC_deterministic_for_GUI(...
                                         weights, params)
% MS 2017.02.23
%
% Function to simulate Wilson-Cowan network activity deterministically at
% the level of excitatory and inhibitory populations
%

% Model weights
W_ee = weights.ee;
W_ie = weights.ie;
W_ei = weights.ei;
W_ii = weights.ii;

% Parameters
T = params.T;
dt = params.dt;
t = (T(1):dt:T(2))';
P = params.P;
Q = params.Q;
alpha_E = params.alpha_E;
alpha_I = params.alpha_I;
theta_E = params.theta_E;
theta_I = params.theta_I;
inact_E = params.inact_E;
inact_I = params.inact_I;

disp(['Weights = ', num2str([W_ee, W_ie, W_ei, W_ii])])

% Initialize vectors
E = zeros(length(t),1); % Mean excitatory membrane potential
I = zeros(length(t),1); % Mean inhibitory membrane potential
E(1) = 0.0;
I(1) = 0.0;
J_E = zeros(length(t),1); % Total synaptic input to excitatory population
J_I = zeros(length(t),1); % Total synaptic input to inhibitory population
F_E = zeros(length(t),1); % Mean firing rate of excitatory population
F_I = zeros(length(t),1); % Mean firing rate of inhibitory population


%%% DRIVER %%%

for i = 1:(length(t)-1)
    J_E(i) = W_ee*E(i) - W_ei*I(i) + P;
    J_I(i) = W_ie*E(i) - W_ii*I(i) + Q;
    F_E(i) = sigmoid_firing_rate(J_E(i), alpha_E, theta_E);
    F_I(i) = sigmoid_firing_rate(J_I(i), alpha_I, theta_I);
    
    dEdt = -1 * inact_E * E(i) + (1 - E(i)) * F_E(i);
    dIdt = -1 * inact_I * I(i) + (1 - I(i)) * F_I(i);
    E(i+1) = E(i) + dt*dEdt;
    I(i+1) = I(i) + dt*dIdt;
end


%%%  ACTIVITY MODEL FUNCTIONS %%%

    function [output] = sigmoid_firing_rate(X,alpha,theta)
        % Sigmoid firing rate function (Wilson and Cowan, 1972)
        %   Parameters set location and value of maximum slope:
        %   max[S'(X)] = S'(theta) = alpha/4
        output = (1/(1+exp(-alpha*(X-theta)))) - (1/(1+exp(alpha*theta)));
    end

    function [output] = gaussian_firing_rate(Jxk,Xsd,Xtheta)
        % Gaussian firing rate function (Meijer et al., 2015)
        output = exp(-1*((Jxk-Xtheta)/Xsd)^2)-exp(-(-Xtheta/Xsd)^2);
    end
end