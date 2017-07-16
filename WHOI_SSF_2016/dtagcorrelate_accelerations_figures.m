% Generate accelerometer cross-correlation plots from 'savedSignals'
% Nicely all auto-correlation plots (3 axes on each plot) and displays peak
% coefficients across the diagonal.
% MS 2016.08.17

N = size(savedSignals,1);

coef = zeros(N,N,3);
lags = zeros(N,N,3);

% Set up figure and axes
close all
f = figure('Position',[50 50 600 500]);
ax = gobjects(N,N);
count = 1;
for row = 1:N
    for col = 1:N
        ax(row,col) = subplot(N,N,count);
        hold(ax(row,col),'on');
        count = count + 1;
    end
end

% Populate axes with cross-correlation plots
colors = ['m','g','b'];
handles = gobjects(N,N,3);
for row = 1:N
    for col = 1:N
        sigRow = savedSignals{row,1}.A_dec;
        sigCol = savedSignals{col,1}.A_dec;
        
        for i = 1:3 % Iterate through X, Y, Z axes
            % Normalize signals by mean-subtracting and dividing by the
            % signals' standard deviation. Normalize for signal length by
            % dividing by 1/sqrt(n). Then calculate cross-correlation of
            % normalized signals
            
            tmp_sigRow = sigRow(:,i);
            tmp_sigCol = sigCol(:,i);
            
            % Normalize signal 1
            tmp_sigRow = (tmp_sigRow - mean(tmp_sigRow))/std(tmp_sigRow);
            % Normalize signal 2
            tmp_sigCol = (tmp_sigCol - mean(tmp_sigCol))/std(tmp_sigCol);
            % Calculate cross-correlation, take abs(), normalize by lengths
            [tmp_xcor, tmp_lag] = xcorr(tmp_sigRow,tmp_sigCol);
            tmp_xcor = abs(tmp_xcor);
            tmp_xcor = tmp_xcor / (sqrt(length(tmp_sigRow))*sqrt(length(tmp_sigCol)));
            
            handles(row,col,i) = plot(ax(row,col),tmp_lag,tmp_xcor,...
                'color',colors(i));
            set(ax(row,col),'xlim',[-1000 1000],'ylim',[0 1],...
                'xtick',[0],'ytick',[0 1])
            
            % Determine the normalized peak cross-correlation
            % coefficient and the offset
            [tmp_coef, tmp_idx] = max(tmp_xcor);
            tmp_offset = tmp_lag(tmp_idx);
            
            coef(row,col,i) = tmp_coef;
            lags(row,col,i) = tmp_offset;
        end
            
    end
end


% Delete repeated cross-correlation plots and visualize coefficients
for row = 1:N
    for col= 1:(row-1)
        set(ax(row,col),'Visible','off')
        delete(handles(row,col,:));
        
        coefStr = {['Coef_{ X }: ',num2str(coef(col,row,1),'%2.2f')];
            ['Coef_{ Y }: ',num2str(coef(col,row,2),'%2.2f')];
            ['Coef_{ Z }: ',num2str(coef(col,row,3),'%2.2f')]};
        
        t = text(0,0.5,coefStr,'Parent',ax(row,col),...
            'Interpreter','TEX','FontWeight','bold',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
        
        if round(min([coef(col,row,1),coef(col,row,2),coef(col,row,3)]),2) >= 0.25
           set(t,'Color','r','FontSize',12)
        end
    end
end

for i = 1:N
    t(1) = title(ax(1,i),['Call ',num2str(i)]);
    set(ax,'YAxislocation','right')
    t(2) = ylabel(ax(i,N),['   Call ',num2str(i)]);
    set(t,'Rotation',0,'FontSize',10,'FontWeight','bold');
end

% % % % Modify auto-correlation plots
% % % for i = 1:N
% % %     set(ax(i,i),'Visible','off')
% % %     delete(handles(i,i,:));
% % % end

set(f,'PaperPositionMode','auto')
%print(f,'TEMP.png','-dpng','-r900')

