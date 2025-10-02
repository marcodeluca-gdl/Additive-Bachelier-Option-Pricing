function plotBachelierParameters(eta, kappa, eta_t, kappa_t, ttm)
% 2 subplots
figure()

% Eta
subplot(2,1,1)
plot(ttm(2:end), eta_t, 'gs--', 'LineWidth', 1.2, 'MarkerSize', 5); hold on;
yline(eta, 'b-', 'LineWidth', 1.2);
xlabel('$t$ [years]', 'Interpreter', 'latex')
ylabel('$\eta$', 'Interpreter', 'latex')
legend({'$\eta_t$', '$\eta$'}, 'Interpreter', 'latex', 'Location', 'northeast')
grid on

% Kappa
subplot(2,1,2)
plot(ttm(2:end), kappa_t, 'gs--', 'LineWidth', 1.2, 'MarkerSize', 5); hold on;
yline(kappa, 'b-', 'LineWidth', 1.2);
xlabel('$t$ [years]', 'Interpreter', 'latex')
ylabel('$k$', 'Interpreter', 'latex')
legend({'$k_t$', '$k$'}, 'Interpreter', 'latex', 'Location', 'northeast')
grid on
end