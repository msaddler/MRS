% Cross-correlation coefficient linear relationships

close all; figure; hold all

%%% Z
subplot(1,3,1)
a = out(:,3);
b = mean(out(:,[1 2]),2);
idx = find(b >= 10);
% a = a(idx);
% b = b(idx);

plot(a,b,'o')
xlim([0 70]); ylim([0 70])

p = polyfit(a,b,1);
fit = p(1)*a + p(2);
hold all; plot(a,fit)

resid = b - fit;
SSresid = sum(resid.^2);
SStotal = (length(b)-1)*var(b);
zrsq = 1 - SSresid/SStotal

%%% Y
subplot(1,3,2)
a = out(:,2);
b = mean(out(:,[1 3]),2);
idx = find(b >= 10);
% a = a(idx);
% b = b(idx);

plot(a,b,'o')
xlim([0 70]); ylim([0 70])

p = polyfit(a,b,1);
fit = p(1)*a + p(2);
hold all; plot(a,fit)

resid = b - fit;
SSresid = sum(resid.^2);
SStotal = (length(b)-1)*var(b);
yrsq = 1 - SSresid/SStotal

%%% X
subplot(1,3,3)
a = out(:,1);
b = mean(out(:,[2 3]),2);
idx = find(b >= 10);
% a = a(idx);
% b = b(idx);

plot(a,b,'o')
xlim([0 70]); ylim([0 70])

p = polyfit(a,b,1);
fit = p(1)*a + p(2);
hold all; plot(a,fit)

resid = b - fit;
SSresid = sum(resid.^2);
SStotal = (length(b)-1)*var(b);
xrsq = 1 - SSresid/SStotal