function plot_zerorates(dates,zero_rates)
figure;
plot(dates(2:end), zero_rates, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
grid on;

title('Zero Rates over Time', 'FontSize', 14);
xlabel('Date', 'FontSize', 12);
ylabel('Zero Rate', 'FontSize', 12);

datetick('x', 'yyyy', 'keeplimits');
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
set(gcf, 'Color', 'w');
box on;