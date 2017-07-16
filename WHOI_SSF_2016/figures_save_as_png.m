% Script for saving high resolution .png images for poster
% MS 2016.08.03

f = gcf;
set(f,'PaperPositionMode','auto')
print(f,'bm15_054a_D_calls_on_depth_profile','-dpng','-r900')
clear f