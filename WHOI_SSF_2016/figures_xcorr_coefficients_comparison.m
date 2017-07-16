

% Load-in signalParams
clear; close all

ox = []; oy = []; oz = []; aud_pk2pk = [];
load('G:\Mark\data_bm\dtag_signals\signalparams\bm15_048a_DS_signalparams_v6_dtagcorrelate')
ox = [ox;cell2mat(signalParams.ox_coef)];
oy = [oy;cell2mat(signalParams.oy_coef)];
oz = [oz;cell2mat(signalParams.oz_coef)];
aud_pk2pk = [aud_pk2pk; cell2mat(signalParams.aud_P2P)];
load('G:\Mark\data_bm\dtag_signals\signalparams\bm15_054a_DS_signalparams_v6_dtagcorrelate')
ox = [ox;cell2mat(signalParams.ox_coef)];
oy = [oy;cell2mat(signalParams.oy_coef)];
oz = [oz;cell2mat(signalParams.oz_coef)];
aud_pk2pk = [aud_pk2pk; cell2mat(signalParams.aud_P2P)];

omax = max([ox,oy,oz],[],2);

focalIdx = [];
% Circle focal D calls (only works with one deployment at a time
% focalIdx = [17 35 62 70 71 72 73 74 75 80];

fig = figure('Position',[75 100 1400 700]);
ax = gobjects(1,4);
for i = 1:length(ax)
    x = 0.06 +(i-1)*0.24;
    y = 0.55;
    w = 0.18;
    h = 0.40;
    ax(i) = axes('Parent',fig,'Position',[x,y,w,h]);
end
%ax(end) = axes('Parent',pan,'Position',[x,0.10,w,0.35]);

xline = [0, .1]; yline = [0.25, 0.25];
for i = 1:3
    plot(ax(i),xline,yline,'k--');
    hold(ax(i),'on')
    xlabel(ax(i),'Audio peak-to-peak amp. (mV)')
end
set(ax(1:3),'xlim',[0 0.1],'ylim',[0 1],'ytick',[0:0.25:1])

plot(ax(4),[0.25 0.25],[0 20],'k--'); hold(ax(4),'on')
histogram(ax(4),omax,20)
set(ax(4),'xlim',[0 1],'xtick',[0:0.25:1],'ylim',[0 16])
xlabel(ax(4),'Max Coefficient')
ylabel(ax(4),'Count')

handles = gobjects(1,3);
handles(1) = plot(ax(1),aud_pk2pk,ox,'k.');
handles(2) = plot(ax(2),aud_pk2pk,oy,'k.');
handles(3) = plot(ax(3),aud_pk2pk,oz,'k.');
set(handles,'MarkerSize',12)
ylabel(ax(1),'Cross-correlation coef. (accel. X-axis)')
ylabel(ax(2),'Cross-correlation coef. (accel. Y-axis)')
ylabel(ax(3),'Cross-correlation coef. (accel. Z-axis)')


if ~isempty(focalIdx)
    h(1) = plot(ax(1),aud_pk2pk(focalIdx),ox(focalIdx),'ko');
    h(2) = plot(ax(2),aud_pk2pk(focalIdx),oy(focalIdx),'ko');
    h(3) = plot(ax(3),aud_pk2pk(focalIdx),oz(focalIdx),'ko');
end

set(fig,'PaperPositionMode','auto')
%print(fig,'TEMP.png','-dpng','-r900')
