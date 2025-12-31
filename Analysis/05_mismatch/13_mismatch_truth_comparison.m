%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Quick Look: True Mismatch Results
%% Comparing manager's decision vs reality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

% Load results
data = readtable('Results/TrueMismatch_Results/true_mismatch_summary.csv');

fprintf('\n========================================\n');
fprintf('TRUE MISMATCH ANALYSIS\n');
fprintf('Manager thinks STABLE, Reality is DECREASE\n');
fprintf('========================================\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Decision Comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('1. DECISION COMPARISON:\n');
fprintf('---------------------------\n');
fprintf('%-15s %-10s %-12s %-12s %-12s\n', ...
    'Month', 'Scenario', 'Manager', 'Correct', 'Error');
fprintf('---------------------------\n');
for i = 1:height(data)
    fprintf('%-15s %-10s %-12d %-12d %-12d\n', ...
        data.month{i}, data.scenario{i}, ...
        data.manager_decision(i), data.correct_decision(i), ...
        data.decision_error(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Population Consequences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n2. FEMALE POPULATION CONSEQUENCES:\n');
fprintf('---------------------------\n');
fprintf('%-15s %-10s %-12s %-12s %-12s %-12s\n', ...
    'Month', 'Scenario', 'Assumed', 'Reality', 'Correct', 'Loss');
fprintf('---------------------------\n');
for i = 1:height(data)
    fprintf('%-15s %-10s %-12.3f %-12.3f %-12.3f %-12.3f\n', ...
        data.month{i}, data.scenario{i}, ...
        data.manager_female_assumed(i), data.reality_female(i), ...
        data.correct_female(i), data.female_loss(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Utility Loss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n3. UTILITY CONSEQUENCES:\n');
fprintf('---------------------------\n');
fprintf('%-15s %-10s %-15s %-15s %-15s %-15s\n', ...
    'Month', 'Scenario', 'Assumed', 'Reality', 'Correct', 'Loss');
fprintf('---------------------------\n');
for i = 1:height(data)
    fprintf('%-15s %-10s %-15.3f %-15.3f %-15.3f %-15.3f\n', ...
        data.month{i}, data.scenario{i}, ...
        data.manager_utility_assumed(i), data.reality_utility(i), ...
        data.correct_utility(i), data.utility_loss(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Summary Statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n4. SUMMARY STATISTICS:\n');
fprintf('---------------------------\n');

% By scenario
scenarios = unique(data.scenario);
fprintf('\nBy Scenario:\n');
for s = 1:length(scenarios)
    scenario_data = data(strcmp(data.scenario, scenarios{s}), :);
    fprintf('  Scenario %s:\n', scenarios{s});
    fprintf('    Mean decision error: %.2f weeks\n', mean(scenario_data.decision_error));
    fprintf('    Mean female loss: %.3f\n', mean(scenario_data.female_loss));
    fprintf('    Mean utility loss: %.3f\n', mean(scenario_data.utility_loss));
    fprintf('    Max female loss: %.3f\n', max(scenario_data.female_loss));
end

% By month
fprintf('\nBy Decision Timing:\n');
months = unique(data.month, 'stable');
for m = 1:length(months)
    month_data = data(strcmp(data.month, months{m}), :);
    fprintf('  %s:\n', months{m});
    fprintf('    Mean decision error: %.2f weeks\n', mean(month_data.decision_error));
    fprintf('    Mean female loss: %.3f\n', mean(month_data.female_loss));
    fprintf('    Mean utility loss: %.3f\n', mean(month_data.utility_loss));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Key Findings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n========================================\n');
fprintf('KEY FINDINGS:\n');
fprintf('========================================\n\n');

% Worst case scenarios
[max_decision_error, idx] = max(abs(data.decision_error));
fprintf('Largest decision error: %d weeks\n', data.decision_error(idx));
fprintf('  Occurred in: %s, Scenario %s\n', data.month{idx}, data.scenario{idx});
fprintf('  Manager chose: %d weeks, Should have: %d weeks\n\n', ...
    data.manager_decision(idx), data.correct_decision(idx));

[max_female_loss, idx] = max(abs(data.female_loss));
fprintf('Largest female density loss: %.3f\n', data.female_loss(idx));
fprintf('  Occurred in: %s, Scenario %s\n', data.month{idx}, data.scenario{idx});
fprintf('  Manager assumed: %.3f, Reality: %.3f, Should be: %.3f\n\n', ...
    data.manager_female_assumed(idx), data.reality_female(idx), data.correct_female(idx));

[max_utility_loss, idx] = max(data.utility_loss);
fprintf('Largest utility loss: %.3f\n', data.utility_loss(idx));
fprintf('  Occurred in: %s, Scenario %s\n', data.month{idx}, data.scenario{idx});

% Risk assessment
fprintf('\n\nRISK ASSESSMENT:\n');
fprintf('---------------------------\n');
at_risk_count = sum(data.reality_female < 1.0);
if at_risk_count > 0
    fprintf('WARNING: %d scenarios result in female density < 1.0\n', at_risk_count);
    at_risk = data(data.reality_female < 1.0, :);
    for i = 1:height(at_risk)
        fprintf('  %s, Scenario %s: Female = %.3f\n', ...
            at_risk.month{i}, at_risk.scenario{i}, at_risk.reality_female(i));
    end
else
    fprintf('No scenarios result in critical female density (< 1.0)\n');
end

% Percent difference in outcomes
fprintf('\n\nPERCENT IMPACT OF MISMATCH:\n');
fprintf('---------------------------\n');
avg_female_loss_pct = mean(data.female_loss ./ data.correct_female * 100);
avg_utility_loss_pct = mean(data.utility_loss ./ data.correct_utility * 100);
fprintf('Average female density loss: %.1f%%\n', avg_female_loss_pct);
fprintf('Average utility loss: %.1f%%\n', avg_utility_loss_pct);

fprintf('\n========================================\n');