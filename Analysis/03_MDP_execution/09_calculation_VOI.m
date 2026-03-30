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
% September: Both recruitment and mast known 
%            Recruitment = observed from brood surveys
%            Mast = State Variable
%
% VALUE OF INFORMATION METRICS:
% ------------------------------
% Recruitment VOI: ΔEV(Jan→Apr) = EV(April) - EV(January)
%                  Change in expected population state (M, J, F) from
%                  incorporating weather-based recruitment forecasts
%
% Total VOI:       ΔEV(Jan→Sept) = EV(September) - EV(January)
%                  Change in expected population state (M, J, F) from
%                  resolving both recruitment and mast uncertainties
%                  Averaged across warm/cold weather conditions
%
% Mast VOI:        ΔEV(Apr→Sept) = EV(September) - EV(April)
%                  Change in expected population state (M, J, F) from
%                  knowing actual mast on top of predicted recruitment
%                  Averaged across warm/cold weather conditions
%
% STRUCTURE:
% ----------
% Results loaded from two sources:
%   - Jan/April:   Results/Utility_Results_UpdatedPlots/all_utility_results.mat
%   - September:   Results/Sept_Mast_Results_StateVar_avg/avg_all_utility_results.mat
%
% VOI computed across all combinations of:
%   - Population trend:  Increase, Stable, Decrease
%   - Scenario:          A, B, C
%   - Weather condition: warm, cold (April and September only)
%
% Outputs:
%   - Full VOI table printed to console (trend x scenario x weather)
%   - Mean VOI summary by population trend
%   - Bar plot saved to Results/VOI_comparison.png
%
% Interpretation:
% ---------------
% VOI is measured as change in expected steady-state abundance of adult
% males (M), juvenile males (J), and females (F) under the optimal policy
% at each decision point. The management goal is to balance harvest
% opportunity (season length) with long-run population sustainability.
%
% Positive (+): Later decision point yields higher expected abundance —
%               additional information enables harvest decisions that
%               better sustain population levels while maintaining
%               hunting opportunity
%
% Zero (0):     Decision timing does not affect expected population state —
%               the optimal harvest decision and resulting abundance are
%               the same regardless of when it is made
%
% Negative (-): Later decision point yields lower expected abundance —
%               decisions made under uncertainty are more precautionary,
%               producing more conservative harvest and higher long-run
%               abundance than decisions made with realized information
%               (e.g., observed high mast in September permits more
%               aggressive harvest, reducing future abundance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load all three result files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load baseline (average mast = 12.998)
load('Results/Sept_Utility_Results_CurrentFuture/avg_CF_all_utility_results.mat');
all_results_avg = all_results;
clear all_results;

% Load low mast (mast = 5)
load('Results/Sept_Utility_Results_CurrentFuture/low_CF_all_utility_results.mat');
all_results_low = all_results;
clear all_results;

% Load high mast (mast = 20)
load('Results/Sept_Utility_Results_CurrentFuture/all_utility_results.mat');
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

% Check weather column for September
if ismember('weather', T_all.Properties.VariableNames)
    fprintf('Unique weather values: %s\n', strjoin(unique(T_all.weather), ', '));
end
fprintf('Unique model variations: %s\n\n', strjoin(unique(T_all.model_var(~strcmp(T_all.model_var, ''))), ', '));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate VOI for each combination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trends = unique(T_all.trend);
scenarios = unique(T_all.scenario);
april_vars = {'Low', 'High'};
weathers = {'warm', 'cold'};  % Add weather dimension

results_voi = [];
manuscript_table_rows = [];  % For simplified manuscript table

for t = 1:length(trends)
    for s = 1:length(scenarios)
        for av = 1:length(april_vars)
            for wth = 1:length(weathers)  % Add weather loop
                trend = trends{t};
                scenario = scenarios{s};
                april_var = april_vars{av};
                weather = weathers{wth};
                
                fprintf('\n========================================\n');
                fprintf('%s - Scenario %s - April Var: %s - Weather: %s\n', ...
                    trend, scenario, april_var, weather);
                fprintf('========================================\n');
                
                % Get January utility (no weather dimension)
                idx_jan = strcmp(T_all.month, 'January') & ...
                          strcmp(T_all.trend, trend) & ...
                          strcmp(T_all.scenario, scenario) & ...
                          strcmp(T_all.mast_condition, 'average');
                
                % Get April utility (specific weather AND variation)
                idx_apr = strcmp(T_all.month, 'April') & ...
                          strcmp(T_all.trend, trend) & ...
                          strcmp(T_all.scenario, scenario) & ...
                          strcmp(T_all.mast_condition, 'average') & ...
                          strcmp(T_all.model_var, april_var) & ...
                          strcmp(T_all.weather, weather);
                
                % Get September utilities (specific weather, each mast condition)
                idx_sep_low = strcmp(T_all.month, 'September') & ...
                              strcmp(T_all.trend, trend) & ...
                              strcmp(T_all.scenario, scenario) & ...
                              strcmp(T_all.mast_condition, 'low') & ...
                              strcmp(T_all.weather, weather);
                
                idx_sep_avg = strcmp(T_all.month, 'September') & ...
                              strcmp(T_all.trend, trend) & ...
                              strcmp(T_all.scenario, scenario) & ...
                              strcmp(T_all.mast_condition, 'average') & ...
                              strcmp(T_all.weather, weather);
                
                idx_sep_high = strcmp(T_all.month, 'September') & ...
                               strcmp(T_all.trend, trend) & ...
                               strcmp(T_all.scenario, scenario) & ...
                               strcmp(T_all.mast_condition, 'high') & ...
                               strcmp(T_all.weather, weather);
                
                % Extract utilities
                U_jan_vec = T_all.expected_utility(idx_jan);
                U_apr_vec = T_all.expected_utility(idx_apr);
                U_sep_low_vec = T_all.expected_utility(idx_sep_low);
                U_sep_avg_vec = T_all.expected_utility(idx_sep_avg);
                U_sep_high_vec = T_all.expected_utility(idx_sep_high);
                
                % Debug
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
                    
                    % Recruitment VOI (Jan → Apr)
                    recruitment_VOI = U_apr - U_jan;
                    
                    % Total VOI (Jan → Sep) for each mast condition
                    total_VOI_low = U_sep_low - U_jan;
                    total_VOI_avg = U_sep_avg - U_jan;
                    total_VOI_high = U_sep_high - U_jan;
                    
                    % Mast VOI (Apr → Sep) for each mast condition
                    mast_VOI_low = U_sep_low - U_apr;
                    mast_VOI_avg = U_sep_avg - U_apr;
                    mast_VOI_high = U_sep_high - U_apr;
                    
                    % Store detailed results
                    result = struct();
                    result.trend = trend;
                    result.scenario = scenario;
                    result.april_var = april_var;
                    result.weather = weather;
                    result.U_January = U_jan;
                    result.U_April = U_apr;
                    result.U_Sep_low = U_sep_low;
                    result.U_Sep_avg = U_sep_avg;
                    result.U_Sep_high = U_sep_high;
                    
                    result.recruitment_VOI = recruitment_VOI;
                    result.total_VOI_low = total_VOI_low;
                    result.total_VOI_avg = total_VOI_avg;
                    result.total_VOI_high = total_VOI_high;
                    
                    result.mast_VOI_low = mast_VOI_low;
                    result.mast_VOI_avg = mast_VOI_avg;
                    result.mast_VOI_high = mast_VOI_high;
                    
                    result.percent_recruitment = (recruitment_VOI / abs(U_jan)) * 100;
                    result.percent_total_low = (total_VOI_low / abs(U_jan)) * 100;
                    result.percent_total_avg = (total_VOI_avg / abs(U_jan)) * 100;
                    result.percent_total_high = (total_VOI_high / abs(U_jan)) * 100;
                    
                    results_voi = [results_voi; result];
                    
                    % Display results
                    fprintf('  U(January):              %.2f (baseline)\n', U_jan);
                    fprintf('  U(April - %s var, %s):   %.2f (ΔU = %.2f, %.2f%%)\n', ...
                        april_var, weather, U_apr, recruitment_VOI, result.percent_recruitment);
                    fprintf('  ---\n');
                    fprintf('  U(Sep | low mast, %s):   %.2f (ΔU = %.2f, %.2f%%)\n', ...
                        weather, U_sep_low, total_VOI_low, result.percent_total_low);
                    fprintf('  U(Sep | avg mast, %s):   %.2f (ΔU = %.2f, %.2f%%)\n', ...
                        weather, U_sep_avg, total_VOI_avg, result.percent_total_avg);
                    fprintf('  U(Sep | high mast, %s):  %.2f (ΔU = %.2f, %.2f%%)\n', ...
                        weather, U_sep_high, total_VOI_high, result.percent_total_high);
                    
                    % Create simplified manuscript table row (only warm weather, for now)
                    if strcmp(weather, 'warm')
                        ms_row = struct();
                        ms_row.TREND = upper(trend);
                        ms_row.Scenario = scenario;
                        ms_row.recruitment_VOI = recruitment_VOI;
                        ms_row.total_VOI_low = total_VOI_low;
                        ms_row.total_VOI_avg = total_VOI_avg;
                        ms_row.total_VOI_high = total_VOI_high;
                        ms_row.mast_VOI_low = mast_VOI_low;
                        ms_row.mast_VOI_avg = mast_VOI_avg;
                        ms_row.mast_VOI_high = mast_VOI_high;
                        manuscript_table_rows = [manuscript_table_rows; ms_row];
                    end
                else
                    warning('Missing data for %s - Scenario %s - April Var %s - Weather %s', ...
                        trend, scenario, april_var, weather);
                end
            end  % weather loop
        end  % april_vars
    end  % scenarios
end  % trends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create summary tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output_dir = 'Results/Utility_Results/';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

if ~isempty(results_voi)
    % Full detailed table
    VOI_table_full = struct2table(results_voi);
    
    fprintf('\n\n========================================\n');
    fprintf('FULL VALUE OF INFORMATION SUMMARY\n');
    fprintf('========================================\n');
    disp(VOI_table_full);
    
    writetable(VOI_table_full, fullfile(output_dir, 'VOI_summary_full.csv'));
    fprintf('\nFull table saved to: %s\n', fullfile(output_dir, 'VOI_summary_full.csv'));
end

if ~isempty(manuscript_table_rows)
    % Simplified manuscript table (format matching your original)
    MS_table = struct2table(manuscript_table_rows);
    
    % Reorder columns to match your format
    MS_table = MS_table(:, {'TREND', 'Scenario', 'recruitment_VOI', ...
                            'total_VOI_low', 'total_VOI_avg', 'total_VOI_high', ...
                            'mast_VOI_low', 'mast_VOI_avg', 'mast_VOI_high'});
    
    % Rename columns with proper headers
    MS_table.Properties.VariableNames = {...
        'TREND', 'Scenario', 'recruitment_VOI_Jan_Apr', ...
        'total_VOI_Jan_Sept_LOW', 'total_VOI_Jan_Sept_AVG', 'total_VOI_Jan_Sept_HIGH', ...
        'mast_VOI_Apr_Sept_LOW', 'mast_VOI_Apr_Sept_AVG', 'mast_VOI_Apr_Sept_HIGH'};
    
    fprintf('\n\n========================================\n');
    fprintf('MANUSCRIPT TABLE (Warm Weather Only)\n');
    fprintf('========================================\n');
    disp(MS_table);
    
    writetable(MS_table, fullfile(output_dir, 'VOI_manuscript_table.csv'));
    fprintf('\nManuscript table saved to: %s\n', fullfile(output_dir, 'VOI_manuscript_table.csv'));
    
    % Also create a formatted version with proper column headers for LaTeX/Word
    writetable(MS_table, fullfile(output_dir, 'VOI_manuscript_table_formatted.txt'), ...
               'Delimiter', '\t');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualization - Separated by April Variation AND Weather
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(results_voi)
    for av = 1:length(april_vars)
        for wth = 1:length(weathers)
            april_var = april_vars{av};
            weather = weathers{wth};
            
            % Filter results for this april_var AND weather
            idx_var = strcmp({results_voi.april_var}, april_var) & ...
                      strcmp({results_voi.weather}, weather);
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
            
            %% COMBINED FIGURE
            fig_combined = figure('Position', [100, 100, 1600, 800]);
            sgtitle(sprintf('VOI Analysis - April %s Variation, %s Weather', ...
                           april_var, upper(weather)), ...
                    'FontSize', 16, 'FontWeight', 'bold');
            
            % Subplot 1: Recruitment VOI (Jan → Apr)
            subplot(2,3,1);
            bar([results_subset.recruitment_VOI], 'FaceColor', [0.2 0.4 0.7]);
            set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
            ylabel('Utility Gain', 'FontSize', 12);
            title('Recruitment VOI (ΔU: Jan → Apr)', 'FontSize', 14);
            grid on;
            yline(0, 'r--', 'LineWidth', 1.5);
            
            % Subplot 2: Total VOI - Low Mast (Jan → Sep)
            subplot(2,3,2);
            bar([results_subset.total_VOI_low], 'FaceColor', [0.8 0.3 0.3]);
            set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
            ylabel('Utility Change', 'FontSize', 12);
            title('Total VOI | Low Mast (ΔU: Jan → Sep)', 'FontSize', 14);
            grid on;
            yline(0, 'r--', 'LineWidth', 1.5);
            
            % Subplot 3: Total VOI - Avg Mast (Jan → Sep)
            subplot(2,3,3);
            bar([results_subset.total_VOI_avg], 'FaceColor', [0.4 0.6 0.4]);
            set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
            ylabel('Utility Change', 'FontSize', 12);
            title('Total VOI | Avg Mast (ΔU: Jan → Sep)', 'FontSize', 14);
            grid on;
            yline(0, 'r--', 'LineWidth', 1.5);
            
            % Subplot 4: Total VOI - High Mast (Jan → Sep)
            subplot(2,3,4);
            bar([results_subset.total_VOI_high], 'FaceColor', [0.6 0.4 0.8]);
            set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
            ylabel('Utility Change', 'FontSize', 12);
            title('Total VOI | High Mast (ΔU: Jan → Sep)', 'FontSize', 14);
            grid on;
            yline(0, 'r--', 'LineWidth', 1.5);
            
            % Subplot 5: Mast VOI - All Conditions (Apr → Sep)
            subplot(2,3,5);
            hold on;
            b1 = bar(x - 0.25, [results_subset.mast_VOI_low], 0.25, 'FaceColor', [0.8 0.3 0.3]);
            b2 = bar(x, [results_subset.mast_VOI_avg], 0.25, 'FaceColor', [0.4 0.6 0.4]);
            b3 = bar(x + 0.25, [results_subset.mast_VOI_high], 0.25, 'FaceColor', [0.6 0.4 0.8]);
            set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
            ylabel('Utility Change', 'FontSize', 12);
            title('Mast VOI (ΔU: Apr → Sep)', 'FontSize', 14);
            legend({'Low Mast', 'Avg Mast', 'High Mast'}, 'Location', 'best');
            grid on;
            yline(0, 'r--', 'LineWidth', 1.5);
            hold off;
            
            % Subplot 6: Comparison of Recruitment vs Mast contributions
            subplot(2,3,6);
            hold on;
            bar(x - 0.15, [results_subset.recruitment_VOI], 0.3, 'FaceColor', [0.2 0.4 0.7]);
            bar(x + 0.15, [results_subset.mast_VOI_avg], 0.3, 'FaceColor', [0.4 0.6 0.4]);
            set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45);
            ylabel('Utility Change', 'FontSize', 12);
            title('Recruitment vs Mast VOI (Avg Mast)', 'FontSize', 14);
            legend({'Recruitment (Jan→Apr)', 'Mast (Apr→Sep)'}, 'Location', 'best');
            grid on;
            yline(0, 'r--', 'LineWidth', 1.5);
            hold off;
            
            % Save figure
            saveas(fig_combined, fullfile(output_dir, ...
                sprintf('VOI_all_plots_%s_%s.png', april_var, weather)));
            saveas(fig_combined, fullfile(output_dir, ...
                sprintf('VOI_all_plots_%s_%s.pdf', april_var, weather)));
            fprintf('Combined figure saved for April %s variation, %s weather\n', ...
                april_var, weather);
            close(fig_combined);
        end  % weather
    end  % april_vars
end

fprintf('\n========================================\n');
fprintf('Analysis Complete!\n');
fprintf('========================================\n');