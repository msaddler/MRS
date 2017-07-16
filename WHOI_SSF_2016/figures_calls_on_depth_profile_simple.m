% Script to plot calls on the whale's dive profile
% MS 2016.07.15

% Load in 'savedSignals','AUD','ACC'

T = []; % Initialize depth trace time vector
p = []; % Initialize depth trace vector
for i=1:size(ACC.T,1) % Assemble time and depth vectors from all chips
   T = [T, ACC.T{i}];
   p = [p; ACC.p{i}];
end

if length(T) ~= length(p) % Ensure time and depth vectors are of equal length
   len = min(length(T),length(p));
   T = T(1:len);
   p = p(1:len);
end

T = T(1:100:end); % Downsample time vector to reasonable resolution
p = p(1:100:end); % Downsample depth vector to reasonable resolution

T_sig = []; % Initialize time vector for signals
p_sig = []; % Initialize depth vector for signals
for i=1:size(savedSignals,1)
    T_sig = [T_sig,savedSignals{i,2}.T(1)];
    p_sig = [p_sig,mean(savedSignals{i,2}.p)];
end

% Plot the figure
close all
f=figure;
set(gca,'Ydir','reverse')
hold all
plot(T,p,'k')
plot(T_sig,p_sig,'r.')
plot(T_sig,p_sig,'ro')
xlim([T(1) T(end)])
ylim([-5,max(p)+5])
datetick('x','HH:MM','keepticks','keeplimits')
xlabel('Local Time','FontSize',15)
ylabel('Depth (m)','FontSize',15)
set(gca,'FontSize',15)
