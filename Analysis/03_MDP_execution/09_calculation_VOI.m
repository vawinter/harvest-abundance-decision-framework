%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Value of Information Analysis - Expected Value Function Comparison
%% Compares E_π[V(s)] across January, April (warm/Low, warm/High), September
%% VOI = difference in expected value functions across decision points
%% Each E_π[V] weighted by its own stationary distribution
%% 
%% Following Runge (2021) framework
%% Winter et al. 20XX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load January/April results
load('C:/Users/vaw5154/OneDrive - The Pennsylvania State University/PhD/PSUTurkey/turkey_SDP/Results/Utility_Results_UpdatedPlots_JanApr_testingfun/all_utility_results.mat');
jan_apr_results = all_results;

% Load September results (average mast)
load('C:/Users/vaw5154/OneDrive - The Pennsylvania State University/PhD/PSUTurkey/turkey_SDP/Results/Sept_Utility_Results_low_updated/all_utility_results_low.mat');
sept_results = all_results;

% Combine
all_combined = [jan_apr_results sept_results];

% Build lookup table from expected_utility (already E_π[V] from saved scripts)
n = length(all_combined);
comparison = struct();
for i = 1:n
    comparison(i).run_id           = all_combined(i).run_id;
    comparison(i).month            = all_combined(i).month;
    comparison(i).trend            = all_combined(i).trend;
    comparison(i).scenario         = all_combined(i).scenario;
    comparison(i).weather          = all_combined(i).weather;
    comparison(i).model_var        = all_combined(i).model_var;
    comparison(i).expected_utility = all_combined(i).expected_utility;
end
T = struct2table(comparison);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build VOI Table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trends    = {'Increase', 'Stable', 'Decrease'};
scenarios = {'A', 'B', 'C'};

% Storage
row_trend           = {};
row_scenario        = {};
jan_ev              = [];
apr_low_ev          = [];
apr_high_ev         = [];
sept_ev             = [];
delta_apr_low_jan   = [];
delta_apr_high_jan  = [];
delta_sept_jan      = [];
delta_sept_apr_low  = [];
delta_sept_apr_high = [];

for t = 1:length(trends)
    for s = 1:length(scenarios)
        trend    = trends{t};
        scenario = scenarios{s};

        % January
        jan_mask = strcmp(T.month, 'January') & ...
                   strcmp(T.trend, trend) & ...
                   strcmp(T.scenario, scenario);

        % April - warm, Low precision
        apr_low_mask = strcmp(T.month, 'April') & ...
                       strcmp(T.trend, trend) & ...
                       strcmp(T.scenario, scenario) & ...
                       strcmp(T.weather, 'warm') & ...
                       strcmp(T.model_var, 'Low');

        % April - warm, High precision
        apr_high_mask = strcmp(T.month, 'April') & ...
                        strcmp(T.trend, trend) & ...
                        strcmp(T.scenario, scenario) & ...
                        strcmp(T.weather, 'warm') & ...
                        strcmp(T.model_var, 'High');

        % September - warm
        sept_mask = strcmp(T.month, 'September') & ...
                    strcmp(T.trend, trend) & ...
                    strcmp(T.scenario, scenario) & ...
                    strcmp(T.weather, 'warm');

        % Extract E_π[V] for each decision point
        jan_val      = T.expected_utility(jan_mask);
        apr_low_val  = T.expected_utility(apr_low_mask);
        apr_high_val = T.expected_utility(apr_high_mask);
        sept_val     = T.expected_utility(sept_mask);

        % Skip if any missing
        if isempty(jan_val) || isempty(apr_low_val) || ...
           isempty(apr_high_val) || isempty(sept_val)
            fprintf('Missing data for %s-%s, skipping\n', trend, scenario);
            continue;
        end

        % Store E_π[V] values
        row_trend{end+1}    = trend;
        row_scenario{end+1} = scenario;
        jan_ev(end+1)       = jan_val(1);
        apr_low_ev(end+1)   = apr_low_val(1);
        apr_high_ev(end+1)  = apr_high_val(1);
        sept_ev(end+1)      = sept_val(1);

        % VOI differences: ΔE[V] = E_π[V^later] - E_π[V^earlier]
        delta_apr_low_jan(end+1)   = apr_low_val(1)  - jan_val(1);
        delta_apr_high_jan(end+1)  = apr_high_val(1) - jan_val(1);
        delta_sept_jan(end+1)      = sept_val(1)     - jan_val(1);
        delta_sept_apr_low(end+1)  = sept_val(1)     - apr_low_val(1);
        delta_sept_apr_high(end+1) = sept_val(1)     - apr_high_val(1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Print Table to Console
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nValue of Information: Expected Value Function Comparison\n');
fprintf('E_pi[V] for each decision point and VOI differences\n');
fprintf('April scenarios: warm weather, Low and High recruitment precision\n');
fprintf('September: average mast (omu = 12.998), warm weather\n\n');

fprintf('%-12s %-10s | %-8s %-8s %-8s %-8s | %-12s %-12s %-10s %-14s %-14s\n', ...
    'Trend', 'Scenario', ...
    'Jan', 'Apr-Lo', 'Apr-Hi', 'Sept', ...
    'AprLo-Jan', 'AprHi-Jan', 'Sept-Jan', 'Sept-AprLo', 'Sept-AprHi');
fprintf('%s\n', repmat('-', 1, 125));

for i = 1:length(row_trend)
    fprintf('%-12s %-10s | %-8.2f %-8.2f %-8.2f %-8.2f | %-12.2f %-12.2f %-10.2f %-14.2f %-14.2f\n', ...
        row_trend{i}, row_scenario{i}, ...
        jan_ev(i), apr_low_ev(i), apr_high_ev(i), sept_ev(i), ...
        delta_apr_low_jan(i), delta_apr_high_jan(i), delta_sept_jan(i), ...
        delta_sept_apr_low(i), delta_sept_apr_high(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Summary by Trend (mean across scenarios A-C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\nMean VOI by Population Trend (averaged across scenarios A, B, C)\n');
fprintf('%-12s | %-8s %-8s %-8s %-8s | %-12s %-12s %-10s %-14s %-14s\n', ...
    'Trend', ...
    'Jan', 'Apr-Lo', 'Apr-Hi', 'Sept', ...
    'AprLo-Jan', 'AprHi-Jan', 'Sept-Jan', 'Sept-AprLo', 'Sept-AprHi');
fprintf('%s\n', repmat('-', 1, 125));

for t = 1:length(trends)
    trend = trends{t};
    idx = strcmp(row_trend, trend);
    if any(idx)
        fprintf('%-12s | %-8.2f %-8.2f %-8.2f %-8.2f | %-12.2f %-12.2f %-10.2f %-14.2f %-14.2f\n', ...
            trend, ...
            mean(jan_ev(idx)), ...
            mean(apr_low_ev(idx)), ...
            mean(apr_high_ev(idx)), ...
            mean(sept_ev(idx)), ...
            mean(delta_apr_low_jan(idx)), ...
            mean(delta_apr_high_jan(idx)), ...
            mean(delta_sept_jan(idx)), ...
            mean(delta_sept_apr_low(idx)), ...
            mean(delta_sept_apr_high(idx)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save CSV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results_table = table(...
    row_trend', row_scenario', ...
    jan_ev', apr_low_ev', apr_high_ev', sept_ev', ...
    delta_apr_low_jan', delta_apr_high_jan', delta_sept_jan', ...
    delta_sept_apr_low', delta_sept_apr_high', ...
    'VariableNames', {...
    'Trend', 'Scenario', ...
    'Jan_EV', 'AprLow_EV', 'AprHigh_EV', 'Sept_EV', ...
    'Delta_AprLow_Jan', 'Delta_AprHigh_Jan', 'Delta_Sept_Jan', ...
    'Delta_Sept_AprLow', 'Delta_Sept_AprHigh'});

writetable(results_table, 'Results/VOI_ValueFunction_table.csv');
fprintf('\nTable saved to Results/VOI_ValueFunction_table.csv\n');