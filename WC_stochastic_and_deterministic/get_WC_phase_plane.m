function f = get_WC_phase_plane(weights, params)

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

% VECTORS AND ARRAYS
E = 0:0.01:1;
I = 0:0.01:1;
dE = zeros(length(E), length(I));
dI = zeros(length(E), length(I));

% Populate dE and dI arrays
for itr_E = 1:length(E)
    for itr_I = 1:length(I)
        tmp_E = E(itr_E);
        tmp_I = I(itr_I);
        
        J_E = W_ee * tmp_E - W_ei * tmp_I + P;
        J_I = W_ie * tmp_E - W_ii * tmp_I + Q;
        F_E = sigmoid_firing_rate(J_E, alpha_E, theta_E);
        F_I = sigmoid_firing_rate(J_I, alpha_I, theta_I);
        
        dE(itr_E, itr_I) = -1 * inact_E * tmp_E + (1 - tmp_E) * F_E;
        dI(itr_E, itr_I) = -1 * inact_I * tmp_I + (1 - tmp_I) * F_I;
    end
end

% OUTPUT FIGURE
f = figure();
hold all
contour(I,E,dE',[0 0],'g');
contour(I,E,dI',[0 0],'r');

step=3; % make the mesh grid sparser for display purposes
[Ig,Eg]=meshgrid(I,E);  % this creates the grid for the quiver command
Es=Eg(1:step:length(E),1:step:length(I));
Is=Ig(1:step:length(E),1:step:length(I));
dEs=dE(1:step:length(E),1:step:length(I));
dIs=dI(1:step:length(E),1:step:length(I));
quiver(Es, Is, dEs, dIs,3);

title('Fig. 4 Wilson and Cowan 1972 green - E-null; red - I-null ')
xlabel('E')
ylabel('I')

%%%  ACTIVITY MODEL FUNCTIONS %%%

    function [output] = sigmoid_firing_rate(X,alpha,theta)
        % Sigmoid firing rate function (Wilson and Cowan, 1972)
        %   Parameters set location and value of maximum slope:
        %   max[S'(X)] = S'(theta) = alpha/4
        output = (1/(1+exp(-alpha*(X-theta)))) - (1/(1+exp(alpha*theta)));
    end

end