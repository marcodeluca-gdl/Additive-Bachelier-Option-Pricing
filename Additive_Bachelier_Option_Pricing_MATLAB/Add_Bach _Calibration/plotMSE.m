function plotMSE(ttm, MSE_Add, MSE_fixedttm)
% plotMSE - Plot of MSE for Additive Bachelier and Lévy fixed ttm models

figure()
semilogy(ttm(2:end), MSE_fixedttm, 'g--', 'LineWidth', 1.5); % Lévy fixed ttm (verde)
hold on 
semilogy(ttm(2:end), MSE_Add, 'b--', 'LineWidth', 1.5); % Additive Bachelier (blu)

xlabel('t [years]');
ylabel('MSE');
title('Comparison of MSE across maturities');
xlim([0, 3]);
ylim([1e-7, 1e-3]);
grid on;

legend('Lévy fixed ttm', 'Additive Bachelier', 'Location', 'northeast');
end
