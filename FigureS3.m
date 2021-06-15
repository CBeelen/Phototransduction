%% Plot of the stimuli
figure(1); clf;
yline(0, 'k', 'LineWidth', 1.5)
yline(3*8.5, 'g', 'LineWidth', 1.5)
yline(252, 'r', 'LineWidth', 1.5)
yline(695, 'b', 'LineWidth', 1.5)

%log scaled intensities
line([4.906 4.906], [0 978], 'Color', [0.9290, 0.6940, 0.1250], 'Linewidth', 1.5);
line([8.908 8.908], [0 1506], 'Color', [0.9290, 0.6940, 0.1250], 'Linewidth', 1.5);
line([14.91 14.91], [0 1984], 'Color', [0.9290, 0.6940, 0.1250], 'Linewidth', 1.5);
line([22.91 22.91], [0 2476], 'Color', [0.9290, 0.6940, 0.1250], 'Linewidth', 1.5);
line([32.92 32.92], [0 3709], 'Color', [0.9290, 0.6940, 0.1250], 'Linewidth', 3);

text(0.01,0.33,'9.5 photons*{\mu}m^{-2}','Units','normalized','FontSize',11)
text(0.08,0.28,'1 ms','Units','normalized','FontSize',11)
text(0.15,0.47,'32.1 photons*{\mu}m^{-2}','Units','normalized','FontSize',11)
text(0.22,0.42,'1 ms','Units','normalized','FontSize',11)
text(0.3,0.6,'96.4 photons*{\mu}m^{-2}','Units','normalized','FontSize',11)
text(0.37,0.55,'3 ms','Units','normalized','FontSize',11)
text(0.52,0.72,'299 photons*{\mu}m^{-2}','Units','normalized','FontSize',11)
text(0.59,0.67,'5 ms','Units','normalized','FontSize',11)
text(0.78,0.98,'5104 photons*{\mu}m^{-2}','Units','normalized','FontSize',11)
text(0.85,0.93,'60 ms','Units','normalized','FontSize',11)

xlim([0 38])
ylim([0 4100])
set(gca,'YTick',[], 'FontSize', 15)
legend('BG 0 photons*{\mu}m^{-2}s^{-1}', 'BG 8.5 photons*{\mu}m^{-2}s^{-1}',...
    'BG 252 photons*{\mu}m^{-2}s^{-1}', 'BG 695 photons*{\mu}m^{-2}s^{-1}',...
    'Flashes', 'Location', 'northwest');
xlabel('time (s)')
ylabel('Intensity')
