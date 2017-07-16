% downsweep_call_contours.m
% MS 2016.07.29

close all

fig = figure;
ax = axes('Parent',fig);
for i = 1:length(signalContours.t_audio)
    % Get time and contour vectors for audio and accelerometer signal
    t_aud = signalContours.t_audio{i,1}';
    f_aud = signalContours.f_audio{i,1};
    t_acc = signalContours.t_accel{i,1}';
    f_acc = signalContours.f_accel{i,1};
    
    % Plot audio and accel contours
    h_aud = plot(ax,t_aud,f_aud,'kx'); hold on
    h_acc = plot(ax,t_acc,f_acc,'rx'); hold off
    ylim(ax,[0 150])
    
    %     confirmButton = [];
    %     while strcmp(confirmButton,'Yes') ~= 1
    %         rect = getrect(ax);
    %         tLim = [rect(1), rect(1)+rect(3)];
    %         fLim = [rect(2), rect(2)+rect(4)];
    %
    %         tIdx = find(t_acc >= tLim(1) & t_acc <= tLim(2));
    %         fIdx = find(f_acc >= fLim(1) & f_acc <= fLim(2));
    %
    %         deleteIdx = intersect(tIdx,fIdx);
    %         t_delete = t_acc(deleteIdx);
    %         f_delete = f_acc(deleteIdx);
    %
    %         keepIdx = setdiff(1:length(t_acc),deleteIdx);
    %         t_acc = t_acc(keepIdx);
    %         f_acc = f_acc(keepIdx);
    %
    %         delete(h_acc)
    %         h_acc = plot(ax,t_acc,f_acc,'bx');
    %         plot(ax,t_delete,f_delete,'r.')
    %
    %         confirmButton = questdlg('Finished deleting points?');
    %     end
    
    
    % % %     p = polyfit(t_aud,f_aud,1);
    % % %     f_aud_fit = polyval(p,t_aud);
    
    % Cubic regression for audio contour
    p = polyfit(t_aud,f_aud,3);
    f_aud_fit = p(1)*(t_aud.^3) + p(2)*(t_aud.^2) + p(3)*(t_aud.^1) + p(4);
    f_aud_resid = f_aud - f_aud_fit;
    f_aud_SSresid = sum(f_aud_resid.^2);
    f_aud_SStotal = (length(f_aud)-1)*var(f_aud);
    f_aud_rsq = 1 - f_aud_SSresid/f_aud_SStotal;
    
    % % %     b = mean(f_acc)-p(1)*mean(t_acc);
    % % %     q = [p(1) b];
    
    % % %     q = polyfit(t_acc,f_acc,1);
    % % %     f_acc_fit = polyval(q,t_acc);
    
    % Cubic regression for accel contour
    q = polyfit(t_acc,f_acc,3);
    f_acc_fit =q(1)*(t_acc.^3) + q(2)*(t_acc.^2) + q(3)*(t_acc.^1) + q(4);
    f_acc_resid = f_acc - f_acc_fit;
    f_acc_SSresid = sum(f_acc_resid.^2);
    f_acc_SStotal = (length(f_acc)-1)*var(f_acc);
    f_acc_rsq = 1 - f_acc_SSresid/f_acc_SStotal;
    
    % Plot figure
    hold on;
    h_fit_aud = plot(t_aud,f_aud_fit,'k');
    h_fit_acc = plot(t_acc,f_acc_fit,'r');
    hold off
    
    % Format legend and axes
    fit_aud_str = ['AUDIO:   ','R^2 = ',num2str(f_aud_rsq),'   Coefficients: ',num2str(p)];
    fit_acc_str = ['ACCEL:   ','R^2 = ',num2str(f_acc_rsq),'   Coefficients: ',num2str(q)];
    legend([h_fit_aud,h_fit_acc],fit_aud_str,fit_acc_str)
    
    titleStr = [signalParams.file{i,1},'   ',signalParams.datestr{i,1}];
    title(ax,titleStr)
    ylim([0 150])
    
    saveStr = [saveDir,'Bandpower_',num2str(maxRatio),'_',savedSignals{i,1}.call,'_',datestr(savedSignals{i,1}.T(1),'yyyy-mm-dd_HH-MM-SS'),'.png'];
    set(f,'PaperPositionMode','auto')
    print(f,saveStr,'-dpng','-r300')
end