% Call correlations graph: vertices := calls, edges := accelerometer signal
% cross-correlation coefficients

close all;

N = size(coef.avg,1);
THRESHOLD = 0.25;

pos = rand(N,2);

% theta = 0:1/N:1; theta = theta(1:N)'; theta = 2*pi*theta;
% pos(:,1) = 100*cos(theta);
% pos(:,2) = 100*sin(theta);

[~,rIdx] = sort(rand(1,N));
pos = pos(rIdx,:);

f = figure('Position', [200 100 800 800]);
ax = axes('Parent',f);
hold(ax,'on');

edges = gobjects(N,N);
for itr1 = 1:N
    for itr2 = 1:N
        if min([coef.X(itr1,itr2),coef.Y(itr1,itr2),coef.Z(itr1,itr2)]) > THRESHOLD
            xPos = [pos(itr1,1), pos(itr2,1)];
            yPos = [pos(itr1,2), pos(itr2,2)];
            edges(itr1,itr2) = plot(xPos,yPos,'r');
            set(edges(itr1,itr2),'LineWidth',coef.Z(itr1,itr2)/10)
        end
    end
end

vertices = gobjects(N,1);
for i = 1:N
   vertices(i) = plot(ax,pos(i,1),pos(i,2),'ko','MarkerFaceColor',[.9 .9 .9]);
   set(vertices(i),'MarkerSize',12)
   text(pos(i,1),pos(i,2),num2str(i))
end

