function plotSignal(signal, windows)
    figure;
    plot(windows, signal, 'LineWidth', 4);
    xlabel('Time (s)');
%     ylabel('Kurtosis');
    ylabel('Log Likelihood');
    grid on;
    xlim([0 7]);
    set(gca, 'FontSize', 55);
%     myLegend = legend('E13', 'E14', 'E15', 'E16');
end
