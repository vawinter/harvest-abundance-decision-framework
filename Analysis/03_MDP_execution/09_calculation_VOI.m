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
%            CVw = 1.116 / 4.96 (baseline variation)
%
% April:     Recruitment predicted from weather models, mast uncertain
%            Recruitment = weather-based prediction
%            Mast = stochastic (mean = 12.998, SD = 3.684)
%            CVw = Low (0.5x baseline) or High (2.0x baseline)
%
% September: Both recruitment and mast known from field observations
%            Recruitment = observed from brood surveys
%            Mast = observed from mast surveys
%            CVw = 0 (known)
%
% VALUE OF INFORMATION METRICS:
% ------------------------------
% Jan → Apr: Value of improved recruitment predictions
%            = U(April) - U(January)
%            = Utility gain from weather-based recruitment forecasts
%            = Calculated separately for Low and High model variation
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
%            = Calculated separately for Low and High model variation
%
% Interpretation:
% ---------------
% Positive (+): Additional information increases utility (gain from waiting)
% Zero (0): No benefit to additional information (doesn't change decisions)
% Negative (-): Additional information decreases utility (loss from waiting)
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
fprintf('Unique mast conditions: %s\n', strjoin(unique(T_all.mast_condition), ', '));
fprintf('Unique model variations: %s\n\n', strjoin(unique(T_all.model_var), ', '));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate VOI for each combination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trends = unique(T_all.trend);
scenarios = unique(T_all.scenario);
april_vars = {'Low', 'High'};

results_voi = [];

for t = 1:length(trends)
    for s = 1:length(scenarios)
        for av = 1:length(april_vars)
            trend = trends{t};
            scenario = scenarios{s};
            april_var = april_vars{av};
            
            fprintf('\n========================================\n');
            fprintf('%s - Scenario %s - April Variation: %s\n', trend, scenario, april_var);
            fprintf('========================================\n');
            
            % Get January utility (uses average mast = 12.998, stochastic)
            idx_jan = strcmp(T_all.month, 'January') & ...
                      strcmp(T_all.trend, trend) & ...
                      strcmp(T_all.scenario, scenario) & ...
                      strcmp(T_all.mast_condition, 'average');
            
            % Get April utility (uses average mast = 12.998, stochastic, specific variation)
            idx_apr = strcmp(T_all.month, 'April') & ...
                      strcmp(T_all.trend, trend) & ...
                      strcmp(T_all.scenario, scenario) & ...
                      strcmp(T_all.mast_condition, 'average') & ...
                      strcmp(T_all.model_var, april_var);
            
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
            
            % Extract utilities
            U_jan_vec = T_all.expected_utility(idx_jan);
            U_apr_vec = T_all.expected_utility(idx_apr);
            U_sep_low_vec = T_all.expected_utility(idx_sep_low);
            U_sep_avg_vec = T_all.expected_utility(idx_sep_avg);
            U_sep_high_vec = T_all.expected_utility(idx_sep_high);
            
            % Debug: Check if we have exactly 1 value for each
            fprintf('  Found: Jan=%d, Apr=%d, SepLow=%d, SepAvg=%d, SepHigh=%d\n', ...
                length(U_jan_vec), length(U_apr_vec), length(U_sep_low_vec), ...
                length(U_sep_avg_vec), length(U_sep_high_vec));
            
            % Take first value if multiple
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
                
                % EVPXI = Additional value of knowing actual conditions (Apr → Sep)
                EVPXI_low = U_sep_low - U_apr;
                EVPXI_avg = U_sep_avg - U_apr;
                EVPXI_high = U_sep_high - U_apr;
                
                % Store results
                result = struct();
                result.trend = trend;
                result.scenario = scenario;
                result.april_var = april_var;
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
                fprintf('  U(April - %s var):       %.4f (EVSI = %.4f, %.2f%%)\n', ...
                    april_var, U_apr, EVSI, result.percent_EVSI);
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
                warning('Missing data for %s - Scenario %s - April Var %s', trend, scenario, april_var);
            end
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
%% Visualization - Separated by April Variation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(results_voi)
    for av = 1:length(april_vars)
        april_var = april_vars{av};
        
        % Filter results for this april_var
        idx_var = strcmp({results_voi.april_var}, april_var);
        results_subset = results_voi(idx_var);
        
        if isempty(results_subset)
            continue;
        end
        
        % Create labels
        labels = cell(length(results_subset), 1);
        for i = 1:length(results_subset)
            trend_short = results_subset(i).trend;
            if length(trend_short) > 3
                trend_short = trend_short(1:3);
            end
            labels{i} = sprintf('%s-%s', trend_short, results_subset(i).scenario);
        end
        
        x = 1:length(results_subset);
        width = 0.25;
        
        %% COMBINED FIGURE
        fig_combined = figure('Position', [100, 100, 1600, 800]);
        sgtitle(sprintf('Value of Information Analysis - April %s Variation', april_var), ...
                'FontSize', 16, 'FontWeight', 'bold');
        
        % Subplot 1: EVSI
        subplot(2,3,1);
        bar([results_subset.EVSI], 'FaceColor', [0.2 0.4 0.7]);
        set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
        ylabel('Utility Gain', 'FontSize', 12);
        title('Value of Recruitment Predictions (Jan → Apr)', 'FontSize', 14);
        grid on;
        yline(0, 'r--', 'LineWidth', 1.5);
        
        % Subplot 2: EVPI by mast
        subplot(2,3,2);
        hold on;
        bar(x-width, [results_subset.EVPI_low], width, 'FaceColor', [0.9 0.6 0.6]);
        bar(x, [results_subset.EVPI_avg], width, 'FaceColor', [0.9 0.9 0.6]);
        bar(x+width, [results_subset.EVPI_high], width, 'FaceColor', [0.6 0.9 0.6]);
        hold off;
        set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
        ylabel('EVPI (Utility Gain)', 'FontSize', 12);
        title('EVPI (Jan → Sep)', 'FontSize', 14);
        legend({'Low Mast', 'Avg Mast', 'High Mast'}, 'Location', 'best');
        grid on;
        yline(0, 'r--', 'LineWidth', 1.5);
        
        % Subplot 3: EVPXI
        subplot(2,3,3);
        hold on;
        bar(x-width, [results_subset.EVPXI_low], width, 'FaceColor', [0.9 0.6 0.6]);
        bar(x, [results_subset.EVPXI_avg], width, 'FaceColor', [0.9 0.9 0.6]);
        bar(x+width, [results_subset.EVPXI_high], width, 'FaceColor', [0.6 0.9 0.6]);
        hold off;
        set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
        ylabel('Incremental Utility Gain', 'FontSize', 12);
        title('Additional Value (Apr → Sep)', 'FontSize', 14);
        legend({'Low Mast', 'Avg Mast', 'High Mast'}, 'Location', 'best');
        grid on;
        yline(0, 'r--', 'LineWidth', 1.5);
        
        % Subplot 4: Utilities by timing
        subplot(2,3,4);
        hold on;
        plot(x, [results_subset.U_January], 'o-', 'LineWidth', 2, 'MarkerSize', 8);
        plot(x, [results_subset.U_April], 's-', 'LineWidth', 2, 'MarkerSize', 8);
        plot(x, [results_subset.U_Sep_avg], 'd-', 'LineWidth', 2, 'MarkerSize', 8);
        hold off;
        set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
        ylabel('Utility', 'FontSize', 12);
        title('Utility by Decision Timing', 'FontSize', 14);
        legend({'January', sprintf('April (%s)', april_var), 'September (avg)'}, 'Location', 'best');
        grid on;
        
        % Subplot 5: September utilities
        subplot(2,3,5);
        hold on;
        bar(x-width, [results_subset.U_Sep_low], width, 'FaceColor', [0.9 0.6 0.6]);
        bar(x, [results_subset.U_Sep_avg], width, 'FaceColor', [0.9 0.9 0.6]);
        bar(x+width, [results_subset.U_Sep_high], width, 'FaceColor', [0.6 0.9 0.6]);
        hold off;
        set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
        ylabel('Utility', 'FontSize', 12);
        title('September Utility by Mast', 'FontSize', 14);
        legend({'Low Mast', 'Avg Mast', 'High Mast'}, 'Location', 'best');
        grid on;
        
        % Subplot 6: Percent improvements
        subplot(2,3,6);
        hold on;
        bar(x-width, [results_subset.percent_EVPI_low], width, 'FaceColor', [0.9 0.6 0.6]);
        bar(x, [results_subset.percent_EVPI_avg], width, 'FaceColor', [0.9 0.9 0.6]);
        bar(x+width, [results_subset.percent_EVPI_high], width, 'FaceColor', [0.6 0.9 0.6]);
        hold off;
        set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
        ylabel('% Improvement', 'FontSize', 12);
        title('EVPI as % of Baseline', 'FontSize', 14);
        legend({'Low Mast', 'Avg Mast', 'High Mast'}, 'Location', 'best');
        grid on;
        yline(0, 'r--', 'LineWidth', 1.5);
        
        % Save combined figure
        saveas(fig_combined, fullfile(output_dir, sprintf('VOI_all_plots_%s.png', april_var)));
        saveas(fig_combined, fullfile(output_dir, sprintf('VOI_all_plots_%s.pdf', april_var)));
        fprintf('Combined figure saved for April %s variation\n', april_var);
        close(fig_combined);
    end
end

fprintf('\n========================================\n');
fprintf('Analysis Complete!\n');
fprintf('========================================\n');