function [coef,lags,output] = dtagcorrelate_accelerations(savedSignals)
% MS 2016.08.17
%
%

N = size(savedSignals,1);

coef.X = zeros(N,N);
coef.Y = zeros(N,N);
coef.Z = zeros(N,N);
lags.X = zeros(N,N);
lags.Y = zeros(N,N);
lags.Z = zeros(N,N);

for itr1 = 1:N % Iterate through all pairs of signals
    for itr2 = 1:N
        signal1 = savedSignals{itr1,1}.A_dec;
        signal2 = savedSignals{itr2,1}.A_dec;
        
        for i = 1:size(signal1,2) % Iterate through X, Y, Z axes
            
            tmp_sigRow = signal1(:,i); % get signal 1 at current axis
            tmp_sigCol = signal2(:,i); % get signal 2 at current axis
            
            % Normalize signal 1
            tmp_sigRow = (tmp_sigRow - mean(tmp_sigRow))/std(tmp_sigRow);
            % Normalize signal 2
            tmp_sigCol = (tmp_sigCol - mean(tmp_sigCol))/std(tmp_sigCol);
            % Calculate cross-correlation, take abs(), normalize by lengths
            [tmp_xcor, tmp_lag] = xcorr(tmp_sigRow,tmp_sigCol);
            tmp_xcor = abs(tmp_xcor);
            tmp_xcor = tmp_xcor / (sqrt(length(tmp_sigRow))*sqrt(length(tmp_sigCol)));
            
            % Determine the normalized peak cross-correlation
            % coefficient and the offset
            [tmp_coef, tmp_idx] = max(tmp_xcor);
            tmp_offset = tmp_lag(tmp_idx);
            
            if i == 1 % Store X-axis cross-correlation coefficent
                coef.X(itr1,itr2) = tmp_coef;
                lags.X(itr1,itr2) = tmp_offset;
            elseif i == 2 % Store Y-axis cross-correlation coefficent
                coef.Y(itr1,itr2) = tmp_coef;
                lags.Y(itr1,itr2) = tmp_offset;
            elseif i == 3 % Store Z-axis cross-correlation coefficent
                coef.Z(itr1,itr2) = tmp_coef;
                lags.Z(itr1,itr2) = tmp_offset;
            end
        end
    end
end

% Prepare a matrix of all inter-call cross-correlations for spreadsheet
output = [];
for row = 1:N
    for col = 1:row-1
        xxcoef = coef.X(row,col);
        yycoef = coef.Y(row,col);
        zzcoef = coef.Z(row,col);
        tmpOut = [row,col,xxcoef,yycoef,zzcoef];
        output = [output;tmpOut];
    end
end

coef.avg = (coef.X + coef.Y + coef.Z)/3;

%     figure('Position',[100 300 1800 600])
%     subplot(1,4,1)
%     imagesc(1:N,1:N,coef.X)
%     colormap('bone')
%
%     subplot(1,4,2)
%     imagesc(1:N,1:N,coef.Y)
%     colormap('bone')
%
%     subplot(1,4,3)
%     imagesc(1:N,1:N,coef.Z)
%     colormap('bone')
%
%     subplot(1,4,4)
%     imagesc(1:N,1:N,coef.avg)
%     colormap('bone')

end