% Scratch-work script for plotting cross-correlation coefficients and lags.
% All of the variables need to be manually constructed by copy-pasting in
% data from excel spreadsheets.

focalFlag = zeros(size(comm));

for i=1:length(focalFlag)
   if strcmp(comm{i,1},'FOCAL')
      focalFlag(i) = 1; 
   end
end

sum(focalFlag);

%%% CROSS CORRELATION COEFFICIENTS
close all
f = figure;
subplot(1,3,1); hold on
plot(coef(:,1),'r.')
plot(focalFlag.*coef(:,1),'ro')
xlabel('calls'); ylabel('xcorr peak (audio , accel X)')

subplot(1,3,2); hold on
plot(coef(:,2),'b.')
plot(focalFlag.*coef(:,2),'bo')
xlabel('calls'); ylabel('xcorr peak (audio , accel Y)')

subplot(1,3,3); hold on
plot(coef(:,3),'k.')
plot(focalFlag.*coef(:,3),'ko')
xlabel('calls'); ylabel('xcorr peak (audio , accel Z)')

suptitle('Cross-correlation coefficients (focal calls are circled)')

set(f,'PaperPositionMode','auto')
print(f,'xcorr coef plot.png','-dpng','-r300')

%%% CROSS CORRELATION DELAYS
close all
f = figure;

hold on;
plot(xydelay,'kx')

for i=1:length(focalFlag)
    if focalFlag(i)
        plot(i,xydelay(i),'ko')
    end
end

title('Cross-correlation delays (position of peaks in # samples)')
ylabel('delay (number of samples)')
xlabel('Calls (focal are circled)')

set(f,'PaperPositionMode','auto')
print(f,'xcorr lags plot.png','-dpng','-r300')