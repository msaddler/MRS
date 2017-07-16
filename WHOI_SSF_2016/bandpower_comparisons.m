% MS 2016.08.08

clear; close all

dir2load = 'G:\Mark\data_bm\dtag_signals\signalparams\DS_signalparams_v4\';
file2load = 'bm15_054a_dtagsignalparams_DS_contours_noise_1s_with_bandpower_analysis';
load([dir2load,file2load]);

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

figure; hold all

plot(BP,BP_avg,'bsq','MarkerFaceColor','b')
plot(BP,BP_max,'ro')
xlabel('Audio bandpower')
ylabel('Accel bandpower')


