function MS_annotate_axes(ax,str)
% Function plots string in the top right corner of provided axes

ypoint = get(ax,'ylim'); ypoint = ypoint(2)-abs(0.02*ypoint(2));
xpoint = get(ax,'xlim'); xpoint = xpoint(2)-abs(0.02*xpoint(2));
text(xpoint,ypoint,str,'Parent',ax,...
    'HorizontalAlignment','right',...
    'VerticalAlignment','top')
end