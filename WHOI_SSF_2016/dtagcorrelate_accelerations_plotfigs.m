function [coef,lags] = dtagcorrelate_accelerations_plotfigs(savedSignals)

N = size(savedSignals,1);

coef.X = zeros(N,N);
coef.Y = zeros(N,N);
coef.Z = zeros(N,N);
lags.X = zeros(N,N);
lags.Y = zeros(N,N);
lags.Z = zeros(N,N);

close all
f(1) = figure('Position',[50 50 600 500]);
f(2) = figure('Position',[50 50 600 500]);
f(3) = figure('Position',[50 50 600 500]);

count = 1;

for itr1 = 1:N % Iterate through all pairs of signals
    for itr2 = 1:N
        signal1 = savedSignals{itr1,1}.A_dec;
        signal2 = savedSignals{itr2,1}.A_dec;
        
        for i = 1:size(signal1,2) % Iterate through X, Y, Z axes
            
            % Calculate cross-correlation and normalize by dividing by
            % the mean and subtracting by the variance
            [tmp_xcor, tmp_lag] = xcorr(signal1(:,i),signal2(:,i));
            tmp_xcor = abs(tmp_xcor);
            tmp_xcor = (tmp_xcor/mean(tmp_xcor)) - var(tmp_xcor);
            
            % Determine the normalized peak cross-correlation
            % coefficient and the offset
            [tmp_coef, tmp_idx] = max(tmp_xcor);
            tmp_offset = tmp_lag(tmp_idx);
            
            % Plot cross-correlations
            if itr1 <= itr2
                if i == 1 % Store X-axis cross-correlation coefficent
                    coef.X(itr1,itr2) = tmp_coef;
                    lags.X(itr1,itr2) = tmp_offset;
                    
                    set(0,'CurrentFigure',f(1))
                    subplot(N,N,count)
                    plot(tmp_lag,tmp_xcor)
                    ylim([0 30]); xlim([-1500 1500])
                    text(-1000,25,num2str(tmp_coef))
                elseif i == 2 % Store Y-axis cross-correlation coefficent
                    coef.Y(itr1,itr2) = tmp_coef;
                    lags.Y(itr1,itr2) = tmp_offset;
                    
                    set(0,'CurrentFigure',f(2))
                    subplot(N,N,count)
                    plot(tmp_lag,tmp_xcor)
                    ylim([0 30]); xlim([-1500 1500])
                    text(-1000,25,num2str(tmp_coef))
                elseif i == 3 % Store Z-axis cross-correlation coefficent
                    coef.Z(itr1,itr2) = tmp_coef;
                    lags.Z(itr1,itr2) = tmp_offset;
                    
                    set(0,'CurrentFigure',f(3))
                    subplot(N,N,count)
                    plot(tmp_lag,tmp_xcor)
                    ylim([0 30]); xlim([-1500 1500])
                    text(-1000,25,num2str(tmp_coef))
                end
            end
        end
        
        count = count + 1;
    end
end

coef.avg = (coef.X + coef.Y + coef.Z)/3;

set(0,'CurrentFigure',f(1))
suptitle('Accel X-axis cross-correlations')
set(0,'CurrentFigure',f(2))
suptitle('Accel Y-axis cross-correlations')
set(0,'CurrentFigure',f(3))
suptitle('Accel Z-axis cross-correlations')


set(f(1),'PaperPositionMode','auto')
%print(f(1),'70-74-80-66 xaxis xcorr','-dpng','-r600')
set(f(2),'PaperPositionMode','auto')
%print(f(2),'70-74-80-66 yaxis xcorr','-dpng','-r600')
set(f(3),'PaperPositionMode','auto')
%print(f(3),'70-74-80-66 zaxis xcorr','-dpng','-r600')

end