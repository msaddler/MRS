% Script to plot calls on the whale's dive profile
% MS 2016.07.15

% EDITED FOR SLIPS AND LUNGES (MS 2017.08.18)

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
f=figure('Position',[50 100 1200 300]);
set(gca,'Ydir','reverse')
%rectangle('Position',[T(1),-150,T(end)-T(1),150.2],'EdgeColor','b','FaceColor','b')
hold all
plot(T,p,'k')

%
% EDITS BELOW
%
T_lunge = lunge_idx_10Hz * (1/10);
T_lunge = T_lunge / (60*60*24);
T_lunge = T_lunge + T(1);
p_lunge = zeros(size(T_lunge));
for i = 1:length(T_lunge)
    tmp_idx = find(T <= T_lunge(i), 1, 'last');
    p_lunge(i) = p(tmp_idx);
end
lungeMark = plot(T_lunge, p_lunge, 'sq','MarkerSize',4,...
    'MarkerEdgeColor','g','MarkerFaceColor','g');
%
% EDITS ABOVE
%


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
legend([noAccMark,accMark,focalMark, lungeMark],...
    'D calls not detected on accelerometer',...
    'D calls detected on accelerometer',...
    'D calls with intra-correlated acc. signals',...
    'Feeding lunge',...
    'Location','southeast')
datetick('x','HH:MM','keepticks','keeplimits')
xlabel('Local Time')
ylabel('Depth (m)')
set(gca,'FontSize',14)
set(gca,'Position',[0.05 0.2 .92 .70])
ylim(gca,[-5 75])
xlim(gca, [T_sig(1) - 0.01, T_sig(end) + 0.02])


%
% EDITS BELOW
%
% T_slip = [17538, 31746] * (1/10); % Slip indexes converted to seconds
% T_slip = T_slip / (60*60*24); % Slip times in days;
% T_slip = T_slip + T(1);
% plot(T_slip, [0, 0], 'gsq');


