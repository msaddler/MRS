function growth_model_WC_plot_realtime(ax_handles, t, synCount, radii, ...
                                       W, excIdx, inhIdx)

mean_W_ee = mean(mean(W(excIdx, excIdx))); % 16
mean_W_ie = mean(mean(W(inhIdx, excIdx))); % 18
mean_W_ei = mean(mean(W(excIdx, inhIdx))); % 12
mean_W_ii = mean(mean(W(inhIdx, inhIdx))); % 3


plot(ax_handles(1), t, synCount.E, 'b.')
plot(ax_handles(1), t, synCount.I, 'r.')
plot(ax_handles(2), t, mean(radii(excIdx)), 'b.')
plot(ax_handles(2), t, mean(radii(inhIdx)), 'r.')
plot(ax_handles(3), t, mean_W_ee, 'b.')
plot(ax_handles(3), t, mean_W_ie, 'bo')
plot(ax_handles(3), t, mean_W_ei, 'r.')
plot(ax_handles(3), t, mean_W_ii, 'ro')

drawnow;

end