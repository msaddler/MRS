% Script for plotting audio and acoustic signal bandlimited power,
% especially for looking at relationship to depth
% MS 2016.08.05

clear; close all

dir2load = 'G:\Mark\data_bm\dtag_signals\signalparams\DS_signalparams_v4\';
file2load = 'bm15_054a_dtagsignalparams_DS_contours_noise_1s_with_bandpower_analysis';
load([dir2load,file2load]);

depth = cell2mat(signalParams.depth);
BP = cell2mat(signalParams.audio_bandpower);
BP1 = cell2mat(signalParams.acc1_bandpower);
BP2 = cell2mat(signalParams.acc2_bandpower);
BP3 = cell2mat(signalParams.acc3_bandpower);
BP_avg = mean([BP1,BP2,BP3],2);
BP_max = max([BP1,BP2,BP3],[],2);

ratio1 = BP1./BP;
ratio2 = BP2./BP;
ratio3 = BP3./BP;
ratio_avg = BP_avg./BP;
ratio_max = BP_max./BP;

figure
subplot(1,2,1); hold all
%plot(depth,ratio1,'rx')%BP1)
%plot(depth,ratio2,'rx')%BP2)
%plot(depth,ratio3,'rx')%BP3)
plot(depth,ratio_max,'ro')
plot(depth,ratio_avg,'bsq','MarkerFaceColor','b')

k=1:length(BP);
%text(depth(k),ratio_avg(k),num2str(k'),'Fontsize',6)

xlabel('Depth')
ylabel('Accel bandlimited power / Audio bandlimited power')
title('Strength of accel. signals (DS calls)')
legend('Max accel. axis','Average across 3 axes')


subplot(1,2,2); hold all
plot(depth,BP,'o')
%text(depth(k),BP(k),num2str(k'),'Fontsize',6)
xlabel('Depth')
ylabel('Bandlimited Power of Acoustic Detections')
title('Strength of acoustic signals (DS calls)')


figure;
subplot(1,2,1); hold all
%plot(depth,ratio1,'rx')%BP1)
%plot(depth,ratio2,'rx')%BP2)
%plot(depth,ratio3,'rx')%BP3)
plot(depth,BP_max,'ro')
plot(depth,BP_avg,'bsq','MarkerFaceColor','b')
legend('Max accel. axis','Average across 3 axes')

k=1:length(BP);
%text(depth(k),BP_avg(k),num2str(k'),'Fontsize',6)

xlabel('Depth')
ylabel('Accel bandlimited power')
title('Strength of accel. signals (DS calls)')

subplot(1,2,2); hold all
plot(depth,BP,'o')
%text(depth(k),BP(k),num2str(k'),'Fontsize',6)
xlabel('Depth')
ylabel('Bandlimited Power of Acoustic Detections')
title('Strength of acoustic signals (DS calls)')


