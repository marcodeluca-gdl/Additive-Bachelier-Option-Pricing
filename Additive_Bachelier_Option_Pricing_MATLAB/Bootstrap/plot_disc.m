function plot_disc(dates,disc)

figure;
plot(dates, disc, '-o', 'LineWidth', 2, 'MarkerSize', 6);
grid on;

title('Discount Factor over Time', 'FontSize', 14);
xlabel('Date', 'FontSize', 12);
ylabel('Discount Factor', 'FontSize', 12);

datetick('x', 'yyyy', 'keeplimits');
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
ylim([0.99955, 1.00001]);

set(gcf, 'Color', 'w');
box on;
