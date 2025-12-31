%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Value of Information Analysis - Decision Timing Framework
% Following framework from: Yokota & Thompson 2004, Runge et al. 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script calculates the expected value of resolving different sources
% of uncertainty to evaluate the relative importance of recruitment versus
% mast information for harvest decisions under different population trends.
%
% DECISION POINTS AND INFORMATION AVAILABLE:
% -------------------------------------------
% January:   Both recruitment and mast uncertain
%            Recruitment = prior year estimate
%            Mast = stochastic (mean = 12.998, SD = 3.684)
%
% April:     Recruitment predicted from weather models, mast uncertain
%            Recruitment = weather-based prediction
%            Mast = stochastic (mean = 12.998, SD = 3.684)
%
% September: Both recruitment and mast known from field observations
%            Recruitment = observed from brood surveys
%            Mast = observed from mast surveys
%
% VALUE OF INFORMATION METRICS:
% ------------------------------
% Jan → Apr: Value of improved recruitment predictions
%            = U(April) - U(January)
%            = Utility gain from weather-based recruitment forecasts
%
% Jan → Sep: Expected Value of Perfect Information (EVPI)
%            = U(September|mast) - U(January)
%            = Value of resolving both uncertainties
%            = Calculated separately for low, average, and high mast
%
% Apr → Sep: Additional value of knowing actual conditions
%            = U(September|mast) - U(April)
%            = Value of knowing actual recruitment and mast vs predicted recruitment
%            = Calculated separately for low, average, and high mast
%
% Interpretation:
% ---------------
% Positive (+): Additional information increases utility (gain from waiting)
% Zero (0): No benefit to additional information (doesn't change decisions)
% Negative (-): Additional information decreases utility (loss from waiting)
%
% The utility difference between January and April quantifies the value of
% improved recruitment predictions, while the difference between January and
% September represents the expected value of perfect information (EVPI) for
% resolving both uncertainties. The incremental utility gain from April to
% September quantifies the additional value of knowing actual recruitment
% and mast conditions compared to relying on predicted recruitment while
% mast remains uncertain.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load all three result files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load baseline (average mast = 12.998)
load('Results/Utility_Results/all_utility_results.mat');
all_results_avg = all_results;
clear all_results;

% Load low mast (mast = 5)
load('Results/Utility_Results/LowMast/all_utility_results_LowMast.mat');
all_results_low = all_results;
clear all_results;

% Load high mast (mast = 20)
load('Results/Utility_Results/HighMast/all_utility_results_HighMast.mat');
all_results_high = all_results;
clear all_results;

fprintf('Loaded results:\n');
fprintf('  Average mast: %d runs\n', length(all_results_avg));
fprintf('  Low mast:     %d runs\n', length(all_results_low));
fprintf('  High mast:    %d runs\n\n', length(all_results_high));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert to tables and add mast condition labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_avg = struct2table(all_results_avg);
T_avg.mast_condition = repmat({'average'}, height(T_avg), 1);
T_avg.mast_value = repmat(12.998, height(T_avg), 1);

T_low = struct2table(all_results_low);
T_low.mast_condition = repmat({'low'}, height(T_low), 1);
T_low.mast_value = repmat(5, height(T_low), 1);

T_high = struct2table(all_results_high);
T_high.mast_condition = repmat({'high'}, height(T_high), 1);
T_high.mast_value = repmat(20, height(T_high), 1);

% Combine all data
T_all = [T_avg; T_low; T_high];

fprintf('Combined table size: %d rows\n', height(T_all));
fprintf('Unique months: %s\n', strjoin(unique(T_all.month), ', '));
fprintf('Unique trends: %s\n', strjoin(unique(T_all.trend), ', '));
fprintf('Unique scenarios: %s\n', strjoin(unique(T_all.scenario), ', '));
fprintf('Unique mast conditions: %s\n\n', strjoin(unique(T_all.mast_condition), ', '));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate VOI for each combination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trends = unique(T_all.trend);
scenarios = unique(T_all.scenario);

results_voi = [];

for t = 1:length(trends)
    for s = 1:length(scenarios)
        trend = trends{t};
        scenario = scenarios{s};
        
        fprintf('\n========================================\n');
        fprintf('%s - Scenario %s\n', trend, scenario);
        fprintf('========================================\n');
        
        % Get January utility (uses average mast = 12.998, stochastic)
        idx_jan = strcmp(T_all.month, 'January') & ...
                  strcmp(T_all.trend, trend) & ...
                  strcmp(T_all.scenario, scenario) & ...
                  strcmp(T_all.mast_condition, 'average');
        
        % Get April utility (uses average mast = 12.998, stochastic)  
        idx_apr = strcmp(T_all.month, 'April') & ...
                  strcmp(T_all.trend, trend) & ...
                  strcmp(T_all.scenario, scenario) & ...
                  strcmp(T_all.mast_condition, 'average');
        
        % Get September utilities for each mast condition (known mast)
        idx_sep_low = strcmp(T_all.month, 'September') & ...
                      strcmp(T_all.trend, trend) & ...
                      strcmp(T_all.scenario, scenario) & ...
                      strcmp(T_all.mast_condition, 'low');
        
        idx_sep_avg = strcmp(T_all.month, 'September') & ...
                      strcmp(T_all.trend, trend) & ...
                      strcmp(T_all.scenario, scenario) & ...
                      strcmp(T_all.mast_condition, 'average');
        
        idx_sep_high = strcmp(T_all.month, 'September') & ...
                       strcmp(T_all.trend, trend) & ...
                       strcmp(T_all.scenario, scenario) & ...
                       strcmp(T_all.mast_condition, 'high');
        
        % Extract utilities - CHECK FOR MULTIPLE VALUES
        U_jan_vec = T_all.expected_utility(idx_jan);
        U_apr_vec = T_all.expected_utility(idx_apr);
        U_sep_low_vec = T_all.expected_utility(idx_sep_low);
        U_sep_avg_vec = T_all.expected_utility(idx_sep_avg);
        U_sep_high_vec = T_all.expected_utility(idx_sep_high);
        
        % Debug: Check if we have exactly 1 value for each
        fprintf('  Found: Jan=%d, Apr=%d, SepLow=%d, SepAvg=%d, SepHigh=%d\n', ...
            length(U_jan_vec), length(U_apr_vec), length(U_sep_low_vec), ...
            length(U_sep_avg_vec), length(U_sep_high_vec));
        
        % Take first value if multiple (shouldn't happen!)
        if ~isempty(U_jan_vec), U_jan = U_jan_vec(1); else, U_jan = NaN; end
        if ~isempty(U_apr_vec), U_apr = U_apr_vec(1); else, U_apr = NaN; end
        if ~isempty(U_sep_low_vec), U_sep_low = U_sep_low_vec(1); else, U_sep_low = NaN; end
        if ~isempty(U_sep_avg_vec), U_sep_avg = U_sep_avg_vec(1); else, U_sep_avg = NaN; end
        if ~isempty(U_sep_high_vec), U_sep_high = U_sep_high_vec(1); else, U_sep_high = NaN; end
        
        % Calculate value of information metrics
        if ~isnan(U_jan) && ~isnan(U_apr) && ~isnan(U_sep_low) && ~isnan(U_sep_avg) && ~isnan(U_sep_high)
            
            % EVSI = Value of weather-based recruitment info (Jan → Apr)
            EVSI = U_apr - U_jan;
            
            % EVPI = Value of perfect info (Jan → Sep)
            EVPI_low = U_sep_low - U_jan;
            EVPI_avg = U_sep_avg - U_jan;
            EVPI_high = U_sep_high - U_jan;
            
            % EVPXI = Additional value of mast info (Apr → Sep)
            EVPXI_low = U_sep_low - U_apr;
            EVPXI_avg = U_sep_avg - U_apr;
            EVPXI_high = U_sep_high - U_apr;
            
            % Store results
            result = struct();
            result.trend = trend;
            result.scenario = scenario;
            result.U_January = U_jan;
            result.U_April = U_apr;
            result.U_Sep_low = U_sep_low;
            result.U_Sep_avg = U_sep_avg;
            result.U_Sep_high = U_sep_high;
            
            result.EVSI = EVSI;
            result.EVPI_low = EVPI_low;
            result.EVPI_avg = EVPI_avg;
            result.EVPI_high = EVPI_high;
            
            result.EVPXI_low = EVPXI_low;
            result.EVPXI_avg = EVPXI_avg;
            result.EVPXI_high = EVPXI_high;
            
            result.percent_EVSI = (EVSI / abs(U_jan)) * 100;
            result.percent_EVPI_low = (EVPI_low / abs(U_jan)) * 100;
            result.percent_EVPI_avg = (EVPI_avg / abs(U_jan)) * 100;
            result.percent_EVPI_high = (EVPI_high / abs(U_jan)) * 100;
            
            results_voi = [results_voi; result];
            
            % Display results
            fprintf('  U(January):              %.4f (baseline)\n', U_jan);
            fprintf('  U(April):                %.4f (EVSI = %.4f, %.2f%%)\n', ...
                U_apr, EVSI, result.percent_EVSI);
            fprintf('  ---\n');
            fprintf('  U(Sep | low mast=5):     %.4f (EVPI = %.4f, %.2f%%)\n', ...
                U_sep_low, EVPI_low, result.percent_EVPI_low);
            fprintf('  U(Sep | avg mast=13):    %.4f (EVPI = %.4f, %.2f%%)\n', ...
                U_sep_avg, EVPI_avg, result.percent_EVPI_avg);
            fprintf('  U(Sep | high mast=20):   %.4f (EVPI = %.4f, %.2f%%)\n', ...
                U_sep_high, EVPI_high, result.percent_EVPI_high);
            fprintf('  ---\n');
            fprintf('  EVSI (Jan→Apr):          %.4f\n', EVSI);
            fprintf('  EVPI range:              %.4f to %.4f\n', ...
                min([EVPI_low, EVPI_avg, EVPI_high]), max([EVPI_low, EVPI_avg, EVPI_high]));
            fprintf('  EVPXI range:             %.4f to %.4f\n', ...
                min([EVPXI_low, EVPXI_avg, EVPXI_high]), max([EVPXI_low, EVPXI_avg, EVPXI_high]));
        else
            warning('Missing data for %s - Scenario %s', trend, scenario);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create summary table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(results_voi)
    VOI_table = struct2table(results_voi);
    
    fprintf('\n\n========================================\n');
    fprintf('VALUE OF INFORMATION SUMMARY\n');
    fprintf('========================================\n');
    disp(VOI_table);
    
    % Save table
    output_dir = 'Results/Utility_Results/';
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    writetable(VOI_table, fullfile(output_dir, 'VOI_summary.csv'));
    fprintf('\nTable saved to: %s\n', fullfile(output_dir, 'VOI_summary.csv'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualization - Save Combined and Individual Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(results_voi)
    % Create labels - FIXED VERSION
    labels = cell(length(results_voi), 1);
    for i = 1:length(results_voi)
        trend_short = results_voi(i).trend;
        if length(trend_short) > 3
            trend_short = trend_short(1:3);
        end
        labels{i} = sprintf('%s-%s', trend_short, results_voi(i).scenario);
    end
    
    x = 1:length(results_voi);
    width = 0.25;
    
    %% COMBINED FIGURE - All 6 plots together
    fig_combined = figure('Position', [100, 100, 1600, 800]);
    
    % Subplot 1: Value of improved recruitment predictions
    subplot(2,3,1);
    bar([results_voi.EVSI]);
    set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
    ylabel('Utility Gain', 'FontSize', 12);
    title('Value of Recruitment Predictions (Jan → Apr)', 'FontSize', 14);
    grid on;
    yline(0, 'r--', 'LineWidth', 1.5);
        
    % Subplot 2: EVPI by mast
    subplot(2,3,2);
    hold on;
    bar(x-width, [results_voi.EVPI_low], width, 'FaceColor', [0.9 0.6 0.6]);
    bar(x, [results_voi.EVPI_avg], width, 'FaceColor', [0.9 0.9 0.6]);
    bar(x+width, [results_voi.EVPI_high], width, 'FaceColor', [0.6 0.9 0.6]);
    hold off;
    set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
    ylabel('EVPI (Utility Gain)', 'FontSize', 12);
    title('EVPI (Jan → Sep)', 'FontSize', 14);
    legend({'Low Mast', 'Avg Mast', 'High Mast'}, 'Location', 'best');
    grid on;
    yline(0, 'r--', 'LineWidth', 1.5);
    
    % Subplot 3: Additional value of knowing actual conditions
    subplot(2,3,3);
    hold on;
    bar(x-width, [results_voi.EVPXI_low], width, 'FaceColor', [0.9 0.6 0.6]);
    bar(x, [results_voi.EVPXI_avg], width, 'FaceColor', [0.9 0.9 0.6]);
    bar(x+width, [results_voi.EVPXI_high], width, 'FaceColor', [0.6 0.9 0.6]);
    hold off;
    set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
    ylabel('Incremental Utility Gain', 'FontSize', 12);
    title('Additional Value of Actual Conditions (Apr → Sep)', 'FontSize', 14);
    legend({'Low Mast', 'Avg Mast', 'High Mast'}, 'Location', 'best');
    grid on;
    yline(0, 'r--', 'LineWidth', 1.5);

    % Subplot 4: Utilities by decision timing
    subplot(2,3,4);
    hold on;
    plot(x, [results_voi.U_January], 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    plot(x, [results_voi.U_April], 's-', 'LineWidth', 2, 'MarkerSize', 8);
    plot(x, [results_voi.U_Sep_avg], 'd-', 'LineWidth', 2, 'MarkerSize', 8);
    hold off;
    set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
    ylabel('Utility', 'FontSize', 12);
    title('Utility by Decision Timing', 'FontSize', 14);
    legend({'January', 'April', 'September (avg)'}, 'Location', 'best');
    grid on;
    
    % Subplot 5: September utilities by mast
    subplot(2,3,5);
    hold on;
    bar(x-width, [results_voi.U_Sep_low], width, 'FaceColor', [0.9 0.6 0.6]);
    bar(x, [results_voi.U_Sep_avg], width, 'FaceColor', [0.9 0.9 0.6]);
    bar(x+width, [results_voi.U_Sep_high], width, 'FaceColor', [0.6 0.9 0.6]);
    hold off;
    set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
    ylabel('Utility', 'FontSize', 12);
    title('September Utility by Mast', 'FontSize', 14);
    legend({'Low Mast', 'Avg Mast', 'High Mast'}, 'Location', 'best');
    grid on;
    
    % Subplot 6: Percent improvements
    subplot(2,3,6);
    hold on;
    bar(x-width, [results_voi.percent_EVPI_low], width, 'FaceColor', [0.9 0.6 0.6]);
    bar(x, [results_voi.percent_EVPI_avg], width, 'FaceColor', [0.9 0.9 0.6]);
    bar(x+width, [results_voi.percent_EVPI_high], width, 'FaceColor', [0.6 0.9 0.6]);
    hold off;
    set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
    ylabel('% Improvement', 'FontSize', 12);
    title('EVPI as % of Baseline', 'FontSize', 14);
    legend({'Low Mast', 'Avg Mast', 'High Mast'}, 'Location', 'best');
    grid on;
    yline(0, 'r--', 'LineWidth', 1.5);
    
    % Save combined figure
    saveas(fig_combined, fullfile(output_dir, 'VOI_all_plots.png'));
    saveas(fig_combined, fullfile(output_dir, 'VOI_all_plots.pdf'));
    fprintf('Combined figure saved: VOI_all_plots.png/pdf\n');
    
    
    %% INDIVIDUAL PLOTS - Save each separately
    
    % Plot 1: Value of recruitment predictions
    fig1 = figure('Position', [100, 100, 800, 600]);
    bar([results_voi.EVSI], 'FaceColor', [0.2 0.4 0.7]);
    set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
    ylabel('Utility Gain', 'FontSize', 14);
    title('Value of Recruitment Predictions (Jan → Apr)', 'FontSize', 16);
    grid on;
    yline(0, 'k--', 'LineWidth', 1.5);
    set(gca, 'FontSize', 12);
    saveas(fig1, fullfile(output_dir, 'VOI_1_recruitment_predictions.png'));
    saveas(fig1, fullfile(output_dir, 'VOI_1_recruitment_predictions.pdf'));
    close(fig1);
    
    % Plot 2: EVPI by mast
    fig2 = figure('Position', [100, 100, 800, 600]);
    hold on;
    bar(x-width, [results_voi.EVPI_low], width, 'FaceColor', [0.9 0.4 0.4]);
    bar(x, [results_voi.EVPI_avg], width, 'FaceColor', [0.9 0.8 0.4]);
    bar(x+width, [results_voi.EVPI_high], width, 'FaceColor', [0.4 0.8 0.4]);
    hold off;
    set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
    ylabel('EVPI (Utility Gain)', 'FontSize', 14);
    title('Expected Value of Perfect Information (Jan → Sep)', 'FontSize', 16);
    legend({'Low Mast (5)', 'Avg Mast (13)', 'High Mast (20)'}, 'Location', 'best', 'FontSize', 12);
    grid on;
    yline(0, 'k--', 'LineWidth', 1.5);
    set(gca, 'FontSize', 12);
    saveas(fig2, fullfile(output_dir, 'VOI_2_EVPI.png'));
    saveas(fig2, fullfile(output_dir, 'VOI_2_EVPI.pdf'));
    close(fig2);
    
    % Plot 3: Additional value of actual conditions
    fig3 = figure('Position', [100, 100, 800, 600]);
    hold on;
    bar(x-width, [results_voi.EVPXI_low], width, 'FaceColor', [0.9 0.4 0.4]);
    bar(x, [results_voi.EVPXI_avg], width, 'FaceColor', [0.9 0.8 0.4]);
    bar(x+width, [results_voi.EVPXI_high], width, 'FaceColor', [0.4 0.8 0.4]);
    hold off;
    set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
    ylabel('Incremental Utility Gain', 'FontSize', 14);
    title('Additional Value of Actual Conditions (Apr → Sep)', 'FontSize', 16);
    legend({'Low Mast (5)', 'Avg Mast (13)', 'High Mast (20)'}, 'Location', 'best', 'FontSize', 12);
    grid on;
    yline(0, 'k--', 'LineWidth', 1.5);
    set(gca, 'FontSize', 12);
    saveas(fig3, fullfile(output_dir, 'VOI_3_additional_value.png'));
    saveas(fig3, fullfile(output_dir, 'VOI_3_additional_value.pdf'));
    close(fig3);
    
    % Plot 4: Utilities by decision timing
    fig4 = figure('Position', [100, 100, 800, 600]);
    hold on;
    plot(x, [results_voi.U_January], 'o-', 'LineWidth', 2.5, 'MarkerSize', 10, ...
         'Color', [0.2 0.4 0.7], 'MarkerFaceColor', [0.2 0.4 0.7]);
    plot(x, [results_voi.U_April], 's-', 'LineWidth', 2.5, 'MarkerSize', 10, ...
         'Color', [0.9 0.5 0.2], 'MarkerFaceColor', [0.9 0.5 0.2]);
    plot(x, [results_voi.U_Sep_avg], 'd-', 'LineWidth', 2.5, 'MarkerSize', 10, ...
         'Color', [0.9 0.8 0.2], 'MarkerFaceColor', [0.9 0.8 0.2]);
    hold off;
    set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
    ylabel('Expected Utility', 'FontSize', 14);
    title('Utility by Decision Timing', 'FontSize', 16);
    legend({'January', 'April', 'September (avg mast)'}, 'Location', 'best', 'FontSize', 12);
    grid on;
    set(gca, 'FontSize', 12);
    saveas(fig4, fullfile(output_dir, 'VOI_4_utility_timing.png'));
    saveas(fig4, fullfile(output_dir, 'VOI_4_utility_timing.pdf'));
    close(fig4);
    
    % Plot 5: September utilities by mast
    fig5 = figure('Position', [100, 100, 800, 600]);
    hold on;
    bar(x-width, [results_voi.U_Sep_low], width, 'FaceColor', [0.9 0.4 0.4]);
    bar(x, [results_voi.U_Sep_avg], width, 'FaceColor', [0.9 0.8 0.4]);
    bar(x+width, [results_voi.U_Sep_high], width, 'FaceColor', [0.4 0.8 0.4]);
    hold off;
    set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
    ylabel('Expected Utility', 'FontSize', 14);
    title('September Utility by Mast Condition', 'FontSize', 16);
    legend({'Low Mast (5)', 'Avg Mast (13)', 'High Mast (20)'}, 'Location', 'best', 'FontSize', 12);
    grid on;
    set(gca, 'FontSize', 12);
    saveas(fig5, fullfile(output_dir, 'VOI_5_september_utility.png'));
    saveas(fig5, fullfile(output_dir, 'VOI_5_september_utility.pdf'));
    close(fig5);
    
    % Plot 6: EVPI as % of baseline
    fig6 = figure('Position', [100, 100, 800, 600]);
    hold on;
    bar(x-width, [results_voi.percent_EVPI_low], width, 'FaceColor', [0.9 0.4 0.4]);
    bar(x, [results_voi.percent_EVPI_avg], width, 'FaceColor', [0.9 0.8 0.4]);
    bar(x+width, [results_voi.percent_EVPI_high], width, 'FaceColor', [0.4 0.8 0.4]);
    hold off;
    set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
    ylabel('Relative EVPI (% of Baseline Utility)', 'FontSize', 14);
    title('EVPI as Percent of January Baseline', 'FontSize', 16);
    legend({'Low Mast (5)', 'Avg Mast (13)', 'High Mast (20)'}, 'Location', 'best', 'FontSize', 12);
    grid on;
    yline(0, 'k--', 'LineWidth', 1.5);
    set(gca, 'FontSize', 12);
    saveas(fig6, fullfile(output_dir, 'VOI_6_percent_improvement.png'));
    saveas(fig6, fullfile(output_dir, 'VOI_6_percent_improvement.pdf'));
    close(fig6);
    
    fprintf('\nAll individual plots saved:\n');
    fprintf('  - VOI_1_recruitment_predictions.png/pdf\n');
    fprintf('  - VOI_2_EVPI.png/pdf\n');
    fprintf('  - VOI_3_additional_value.png/pdf\n');
    fprintf('  - VOI_4_utility_timing.png/pdf\n');
    fprintf('  - VOI_5_september_utility.png/pdf\n');
    fprintf('  - VOI_6_percent_improvement.png/pdf\n');
end

fprintf('\n========================================\n');
fprintf('Analysis Complete!\n');
fprintf('========================================\n');