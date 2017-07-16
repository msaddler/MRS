% Script to plot calls on the whale's dive profile
% MS 2016.07.15

% Load in 'savedSignals','AUD','ACC'

focalIdx = [17 35 62 70 71 72 73 74 75 80];

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

% Distinguish between accelerometer detected and not-detected calls
THRESHOLD = 0.25; % <-- threshold value for variable being used
maxCoefs = []; % Initialize vector for distinguishing variable
for i=1:size(savedSignals,1)
   ox_coef = signalParams.ox_coef{i,1};
   oy_coef = signalParams.oy_coef{i,1};
   oz_coef = signalParams.oz_coef{i,1};
   maxCoefs = [maxCoefs; max([ox_coef, oy_coef, oz_coef])];
end

% Plot the figure
close all
f=figure('Position',[50 100 1800 600]);
set(gca,'Ydir','reverse')
%rectangle('Position',[T(1),-150,T(end)-T(1),150.2],'EdgeColor','b','FaceColor','b')
hold all
plot(T,p,'k')
xlim([T(1) T(end)])
ylim([-5,max(p)+5])
for i=1:length(T_sig)
    if maxCoefs(i) >= THRESHOLD
        accMark = plot(T_sig(i),p_sig(i),'o','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r');
    else
        noAccMark = plot(T_sig(i),p_sig(i),'o','MarkerSize',6,'MarkerEdgeColor','b','MarkerFaceColor','b');
    end
    
    if ~isempty(find(focalIdx == i,1))
        focalMark = plot(T_sig(i),p_sig(i),'o','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
    end
end
legend([noAccMark,accMark,focalMark],...
    'D calls not detected on accelerometer',...
    'D calls detected on accelerometer',...
    'D calls with intra-correlated acc. signals',...
    'Location','southeast')
datetick('x','HH:MM','keepticks','keeplimits')
xlabel('Local Time')
ylabel('Depth (m)')
set(gca,'FontSize',14)
set(gca,'Position',[0.05 0.2 .92 .70])
ylim(gca,[-5 75])

