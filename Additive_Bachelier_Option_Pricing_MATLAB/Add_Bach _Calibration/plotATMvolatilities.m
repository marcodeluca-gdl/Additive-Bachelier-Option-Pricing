function plotATMvolatilities(ttm,sigmaATM)
figure;
plot(ttm, sigmaATM, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
grid on;

title('ATM Volatilities vs Time to Maturity', 'FontSize', 14);
xlabel('Time to Maturity (years)', 'FontSize', 12);
ylabel('ATM Volatility', 'FontSize', 12);

set(gca, 'FontSize', 12, 'LineWidth', 1.2);
set(gcf, 'Color', 'w');
box on;