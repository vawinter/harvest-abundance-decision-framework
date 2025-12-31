%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualizing the Expectation Gap: Manager Assumptions vs Reality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

% Load results
data = readtable('Results/TrueMismatch_Results/true_mismatch_summary.csv');

% Create output directory for plots
output_dir = 'Results/TrueMismatch_Results/Plots';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('Creating expectation gap visualizations...\n');

% Define colors
color_assumed = [0.2, 0.6, 0.9];  % Blue - what manager thinks
color_reality = [0.9, 0.3, 0.3];  % Red - what actually happens
color_correct = [0.3, 0.7, 0.3];  % Green - what should be expected

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. The Big Picture: Expectation Gap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', [100, 100, 1200, 500]);

% Female density comparison
subplot(1,2,1);
x = 1:height(data);
width = 0.25;

bar(x - width, data.manager_female_assumed, width, 'FaceColor', color_assumed, 'DisplayName', 'Manager Assumes (Stable)');
hold on;
bar(x, data.reality_female, width, 'FaceColor', color_reality, 'DisplayName', 'Reality (Decrease)');
bar(x + width, data.correct_female, width, 'FaceColor', color_correct, 'DisplayName', 'Correct Expectation');

% Labels
labels = cell(height(data), 1);
for i = 1:height(data)
    if strcmp(data.weather{i}, '')
        labels{i} = sprintf('%s-%s', data.month{i}, data.scenario{i});
    else
        labels{i} = sprintf('%s(%s)-%s', data.month{i}, data.weather{i}(1), data.scenario{i});
    end
end

set(gca, 'XTick', 1:height(data), 'XTickLabel', labels, 'XTickLabelRotation', 45);
ylabel('Female Density', 'FontSize', 12);
title('Female Population: Manager Expectation vs Reality', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
ylim([0, max(data.manager_female_assumed) * 1.1]);
hold off;

% Utility comparison
subplot(1,2,2);
bar(x - width, data.manager_utility_assumed, width, 'FaceColor', color_assumed, 'DisplayName', 'Manager Assumes');
hold on;
bar(x, data.reality_utility, width, 'FaceColor', color_reality, 'DisplayName', 'Reality');
bar(x + width, data.correct_utility, width, 'FaceColor', color_correct, 'DisplayName', 'Correct');

set(gca, 'XTick', 1:height(data), 'XTickLabel', labels, 'XTickLabelRotation', 45);
ylabel('Expected Utility', 'FontSize', 12);
title('Utility: Manager Expectation vs Reality', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
hold off;

sgtitle('The Expectation Gap: Manager Thinks Stable, Reality is Decreasing', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(output_dir, 'expectation_gap_overview.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Percent Overestimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', [100, 100, 800, 600]);

% Calculate percent overestimation
female_overestimate_pct = (data.manager_female_assumed - data.reality_female) ./ data.reality_female * 100;
utility_overestimate_pct = (data.manager_utility_assumed - data.reality_utility) ./ data.reality_utility * 100;

subplot(2,1,1);
bar(x, female_overestimate_pct, 'FaceColor', [0.8, 0.4, 0.4]);
set(gca, 'XTick', 1:height(data), 'XTickLabel', labels, 'XTickLabelRotation', 45);
ylabel('Overestimation (%)', 'FontSize', 12);
title('Female Density: Percent Overestimation by Manager', 'FontSize', 14, 'FontWeight', 'bold');
yline(0, '--k', 'LineWidth', 1.5);
grid on;

subplot(2,1,2);
bar(x, utility_overestimate_pct, 'FaceColor', [0.4, 0.4, 0.8]);
set(gca, 'XTick', 1:height(data), 'XTickLabel', labels, 'XTickLabelRotation', 45);
ylabel('Overestimation (%)', 'FontSize', 12);
title('Utility: Percent Overestimation by Manager', 'FontSize', 14, 'FontWeight', 'bold');
yline(0, '--k', 'LineWidth', 1.5);
grid on;

sgtitle('Percent Overestimation When Manager Assumes Stable Population', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(output_dir, 'percent_overestimation.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. By Decision Timing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', [100, 100, 1000, 600]);

months = {'January', 'April', 'September'};
scenarios = {'A', 'B', 'C'};

for m_idx = 1:length(months)
    month = months{m_idx};
    month_data = data(strcmp(data.month, month), :);
    
    if isempty(month_data)
        continue;
    end
    
    % Female density
    subplot(2, 3, m_idx);
    bar_data = zeros(height(month_data), 3);
    bar_data(:,1) = month_data.manager_female_assumed;
    bar_data(:,2) = month_data.reality_female;
    bar_data(:,3) = month_data.correct_female;
    
    bar(bar_data);
    set(gca, 'XTickLabel', month_data.scenario);
    xlabel('Scenario', 'FontSize', 10);
    ylabel('Female Density', 'FontSize', 10);
    title(sprintf('%s: Female Density', month), 'FontSize', 12, 'FontWeight', 'bold');
    legend({'Assumed', 'Reality', 'Correct'}, 'Location', 'best', 'FontSize', 8);
    grid on;
    
    % Utility
    subplot(2, 3, m_idx + 3);
    bar_data = zeros(height(month_data), 3);
    bar_data(:,1) = month_data.manager_utility_assumed;
    bar_data(:,2) = month_data.reality_utility;
    bar_data(:,3) = month_data.correct_utility;
    
    bar(bar_data);
    set(gca, 'XTickLabel', month_data.scenario);
    xlabel('Scenario', 'FontSize', 10);
    ylabel('Utility', 'FontSize', 10);
    title(sprintf('%s: Utility', month), 'FontSize', 12, 'FontWeight', 'bold');
    legend({'Assumed', 'Reality', 'Correct'}, 'Location', 'best', 'FontSize', 8);
    grid on;
end

sgtitle('Expectation Gap by Decision Timing', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(output_dir, 'expectation_gap_by_timing.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. The Reality Check: What Manager Thinks vs What Happens
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', [100, 100, 800, 800]);

% Female density scatter
subplot(2,2,1);
plot([0, 3], [0, 3], '--k', 'LineWidth', 2, 'DisplayName', 'Perfect Agreement');
hold on;
scatter(data.manager_female_assumed, data.reality_female, 100, 'filled', ...
    'MarkerFaceColor', color_reality, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

xlabel('Manager Assumes (Stable)', 'FontSize', 11);
ylabel('Reality (Decrease)', 'FontSize', 11);
title('Female Density: Assumption vs Reality', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 9);
grid on;
axis square;
hold off;

% Add text annotations for key points
for i = 1:height(data)
    text(data.manager_female_assumed(i), data.reality_female(i), ...
        sprintf(' %s-%s', data.month{i}(1), data.scenario{i}), ...
        'FontSize', 8);
end

% Female density error
subplot(2,2,2);
female_error = data.manager_female_assumed - data.reality_female;
bar(x, female_error, 'FaceColor', [0.9, 0.5, 0.3]);
set(gca, 'XTick', 1:height(data), 'XTickLabel', labels, 'XTickLabelRotation', 45);
ylabel('Overestimation (density units)', 'FontSize', 11);
title('Female Density Overestimation', 'FontSize', 13, 'FontWeight', 'bold');
yline(mean(female_error), '--r', sprintf('Mean = %.3f', mean(female_error)), 'LineWidth', 2, 'FontSize', 10);
grid on;

% Utility scatter
subplot(2,2,3);
plot([100, 260], [100, 260], '--k', 'LineWidth', 2, 'DisplayName', 'Perfect Agreement');
hold on;
scatter(data.manager_utility_assumed, data.reality_utility, 100, 'filled', ...
    'MarkerFaceColor', color_reality, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

xlabel('Manager Assumes (Stable)', 'FontSize', 11);
ylabel('Reality (Decrease)', 'FontSize', 11);
title('Utility: Assumption vs Reality', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 9);
grid on;
axis square;
hold off;

% Add text annotations
for i = 1:height(data)
    text(data.manager_utility_assumed(i), data.reality_utility(i), ...
        sprintf(' %s-%s', data.month{i}(1), data.scenario{i}), ...
        'FontSize', 8);
end

% Utility error
subplot(2,2,4);
utility_error = data.manager_utility_assumed - data.reality_utility;
bar(x, utility_error, 'FaceColor', [0.5, 0.5, 0.9]);
set(gca, 'XTick', 1:height(data), 'XTickLabel', labels, 'XTickLabelRotation', 45);
ylabel('Overestimation (utility units)', 'FontSize', 11);
title('Utility Overestimation', 'FontSize', 13, 'FontWeight', 'bold');
yline(mean(utility_error), '--r', sprintf('Mean = %.1f', mean(utility_error)), 'LineWidth', 2, 'FontSize', 10);
grid on;

sgtitle('Reality Check: Manager Assumptions Systematically Overestimate', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(output_dir, 'reality_check.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Summary Statistics Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', [100, 100, 1000, 400]);

% By scenario
subplot(1,2,1);
scenario_summary = zeros(3, 3);
for s = 1:length(scenarios)
    scenario_data = data(strcmp(data.scenario, scenarios{s}), :);
    scenario_summary(s, 1) = mean(scenario_data.manager_female_assumed);
    scenario_summary(s, 2) = mean(scenario_data.reality_female);
    scenario_summary(s, 3) = mean(scenario_data.correct_female);
end

bar(categorical(scenarios), scenario_summary);
legend({'Manager Assumes', 'Reality', 'Correct'}, 'Location', 'best', 'FontSize', 10);
ylabel('Mean Female Density', 'FontSize', 11);
xlabel('Scenario', 'FontSize', 11);
title('Female Density by Scenario', 'FontSize', 13, 'FontWeight', 'bold');
grid on;

% By month
subplot(1,2,2);
month_summary = zeros(3, 3);
for m = 1:length(months)
    month_data = data(strcmp(data.month, months{m}), :);
    if ~isempty(month_data)
        month_summary(m, 1) = mean(month_data.manager_female_assumed);
        month_summary(m, 2) = mean(month_data.reality_female);
        month_summary(m, 3) = mean(month_data.correct_female);
    end
end

bar(categorical(months), month_summary);
legend({'Manager Assumes', 'Reality', 'Correct'}, 'Location', 'best', 'FontSize', 10);
ylabel('Mean Female Density', 'FontSize', 11);
xlabel('Decision Timing', 'FontSize', 11);
title('Female Density by Decision Timing', 'FontSize', 13, 'FontWeight', 'bold');
grid on;

sgtitle('Summary: Expectation Gap Across Scenarios and Timing', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(output_dir, 'summary_by_groups.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Key Insight Figure: The 46% Gap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', [100, 100, 600, 500]);

% Calculate mean values
mean_assumed = mean(data.manager_female_assumed);
mean_reality = mean(data.reality_female);
mean_gap = mean_assumed - mean_reality;
percent_gap = (mean_gap / mean_reality) * 100;

% Create bar plot
bar_data = [mean_assumed, mean_reality];
b = bar(categorical({'Manager Assumes', 'Reality'}), bar_data, 'FaceWidth', 0.6);
b.FaceColor = 'flat';
b.CData(1,:) = color_assumed;
b.CData(2,:) = color_reality;

% Add value labels on bars
text(1, mean_assumed, sprintf('%.2f', mean_assumed), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontSize', 14, 'FontWeight', 'bold');
text(2, mean_reality, sprintf('%.2f', mean_reality), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Add annotation showing gap
hold on;
plot([1, 2], [mean_assumed, mean_assumed], '--k', 'LineWidth', 1.5);
plot([2, 2], [mean_reality, mean_assumed], 'k', 'LineWidth', 2);
text(2.1, (mean_assumed + mean_reality)/2, ...
    sprintf('%.0f%% Gap\n(%.2f density units)', percent_gap, mean_gap), ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
hold off;

ylabel('Mean Female Density', 'FontSize', 13);
title({'The Expectation Gap:', 'Managers Overestimate Female Population by 46%'}, ...
    'FontSize', 15, 'FontWeight', 'bold');
ylim([0, max(bar_data) * 1.2]);
grid on;

saveas(gcf, fullfile(output_dir, 'key_insight_46percent_gap.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Print Summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n========================================\n');
fprintf('EXPECTATION GAP VISUALIZATIONS COMPLETE\n');
fprintf('========================================\n\n');

fprintf('Key Finding:\n');
fprintf('  Mean female density (Manager assumes): %.2f\n', mean_assumed);
fprintf('  Mean female density (Reality): %.2f\n', mean_reality);
fprintf('  Absolute gap: %.2f density units\n', mean_gap);
fprintf('  Percent overestimation: %.1f%%\n\n', percent_gap);

fprintf('Utility Gap:\n');
fprintf('  Mean utility (Manager assumes): %.1f\n', mean(data.manager_utility_assumed));
fprintf('  Mean utility (Reality): %.1f\n', mean(data.reality_utility));
fprintf('  Utility overestimation: %.1f\n\n', mean(data.manager_utility_assumed - data.reality_utility));

fprintf('All plots saved to: %s\n', output_dir);
fprintf('========================================\n');