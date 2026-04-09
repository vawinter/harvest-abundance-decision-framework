%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualizing the Mismatch Gap: Assumptions vs Reality
%% January 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

% Load results
data = readtable('C:/Users/vaw5154/OneDrive - The Pennsylvania State University/PhD/Manuscripts/Chapter2_harvest-abundance-decision-framework/Results/TrueMismatch_Results_Updated/true_mismatch_summary_JanAp.csv');

% FILTER FOR WARM WEATHER ONLY (or empty for January)
data = data(strcmp(data.weather, 'warm') | strcmp(data.weather, ''), :);

fprintf('\n========================================\n');
fprintf('FILTERED TO WARM WEATHER SCENARIOS ONLY\n');
fprintf('Total rows after filtering: %d\n', height(data));
fprintf('========================================\n\n');

% Create output directory for plots
output_dir = 'C:/Users/vaw5154/OneDrive - The Pennsylvania State University/PhD/Manuscripts/Chapter2_harvest-abundance-decision-framework/Results/TrueMismatch_Results/Plots';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Table 3: Full Contingency Table with Precision Labels (WARM ONLY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create timing labels that include precision for April
timing_labels = cell(height(data), 1);
for i = 1:height(data)
    if strcmp(data.month{i}, 'January')
        timing_labels{i} = 'January';
    elseif strcmp(data.month{i}, 'September')
        timing_labels{i} = 'September';
    elseif strcmp(data.month{i}, 'April')
        if ~isempty(data.precision{i}) && ~strcmp(data.precision{i}, '')
            if strcmp(data.precision{i}, 'High')
                timing_labels{i} = 'April (High)';
            elseif strcmp(data.precision{i}, 'Low')
                timing_labels{i} = 'April (Low)';
            else
                timing_labels{i} = 'April';
            end
        else
            timing_labels{i} = 'April';
        end
    end
end

T3 = table(timing_labels, data.scenario, data.precision, ...
    data.manager_decision, data.correct_decision, ...
    data.manager_female_assumed, data.reality_female, ...
    data.manager_utility_assumed, data.reality_utility, ...
    'VariableNames', {'DecisionTiming', 'HunterScenario', 'Precision', ...
    'ManagerDecision', 'CorrectDecision', ...
    'ExpectedFemale', 'RealizedFemale', ...
    'ExpectedUtility', 'RealizedUtility'});

fprintf('\nTable 3: Complete Results with Precision (Warm Weather Only)\n\n');
disp(T3);

% Export to CSV
writetable(T3, 'Table3_Misidentification_Full_WarmOnly_WithPrecision.csv');
fprintf('Saved: Table3_Misidentification_Full_WarmOnly_WithPrecision.csv\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check for April precision variants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

april_data = data(strcmp(data.month, 'April'), :);
fprintf('\nApril Precision Check:\n');
fprintf('Total April rows: %d\n', height(april_data));
if height(april_data) > 0
    fprintf('Unique precision values: ');
    disp(unique(april_data.precision));
    
    % Count by precision
    for prec = unique(april_data.precision)'
        count = sum(strcmp(april_data.precision, prec{1}));
        fprintf('  Precision "%s": %d rows\n', prec{1}, count);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Summary Statistics (WARM ONLY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

overall_stats = struct();
overall_stats.Cases = height(data);
overall_stats.ManagerDecision = mean(data.manager_decision);
overall_stats.CorrectDecision = mean(data.correct_decision);
overall_stats.ExpectedFemale = mean(data.manager_female_assumed);
overall_stats.RealizedFemale = mean(data.reality_female);
overall_stats.CorrectFemale = mean(data.correct_female);
overall_stats.ExpectedUtility = mean(data.manager_utility_assumed);
overall_stats.RealizedUtility = mean(data.reality_utility);
overall_stats.UtilityLoss_pct = 100 * (1 - overall_stats.RealizedUtility / overall_stats.ExpectedUtility);
overall_stats.ExtinctionCases = sum(data.reality_female < 0.1);

fprintf('\nSummary Statistics (Warm Weather Only):\n');
fprintf('---------------------------\n');
fprintf('Total cases analyzed: %d\n', overall_stats.Cases);
fprintf('Manager decision: %.1f weeks (consistent)\n', overall_stats.ManagerDecision);
fprintf('Correct decision: %.2f weeks (average)\n', overall_stats.CorrectDecision);
fprintf('Expected female density: %.3f\n', overall_stats.ExpectedFemale);
fprintf('Realized female density: %.3f\n', overall_stats.RealizedFemale);
fprintf('Expected utility: %.1f\n', overall_stats.ExpectedUtility);
fprintf('Realized utility: %.1f\n', overall_stats.RealizedUtility);
fprintf('Utility loss: %.1f%%\n', overall_stats.UtilityLoss_pct);
fprintf('Near-extinction cases (F < 0.1): %d/%d\n', overall_stats.ExtinctionCases, overall_stats.Cases);
fprintf('---------------------------\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Manuscript-Ready Formatted Table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n========================================\n');
fprintf('MANUSCRIPT TABLE (WARM WEATHER WITH PRECISION)\n');
fprintf('========================================\n\n');

fprintf('Table X. Management outcomes when declining populations are misidentified as stable\n');
fprintf('(warm weather/high recruitment scenarios).\n\n');

fprintf('%-15s %-10s %-10s %-12s %-12s %-15s %-15s %-15s %-15s\n', ...
    'Decision', 'Hunter', 'Precision', 'Manager', 'Correct', 'Expected', 'Realized', 'Expected', 'Realized');
fprintf('%-15s %-10s %-10s %-12s %-12s %-15s %-15s %-15s %-15s\n', ...
    'Timing', 'Scenario', '', 'Decision', 'Decision', 'Female', 'Female', 'Utility', 'Utility');
fprintf(repmat('-', 1, 130));
fprintf('\n');

for i = 1:height(data)
    prec_str = '';
    if ~isempty(data.precision{i}) && ~strcmp(data.precision{i}, '')
        prec_str = data.precision{i};
    end
    
    fprintf('%-15s %-10s %-10s %-12d %-12.2f %-15.3f %-15.3f %-15.1f %-15.1f\n', ...
        timing_labels{i}, data.scenario{i}, prec_str, ...
        data.manager_decision(i), data.correct_decision(i), ...
        data.manager_female_assumed(i), data.reality_female(i), ...
        data.manager_utility_assumed(i), data.reality_utility(i));
end

fprintf(repmat('-', 1, 130));
fprintf('\n\n');

fprintf('Note: Manager consistently selected 2-week seasons based on stable-population\n');
fprintf('assumptions. Results shown for warm weather conditions (high recruitment predictions).\n');
if any(~strcmp(data.precision, ''))
    fprintf('April results include High and Low precision model variants (different levels of\n');
    fprintf('uncertainty around recruitment estimates).\n');
end
fprintf('Correct decision averages %.2f weeks (primarily closure).\n', overall_stats.CorrectDecision);
fprintf('All %d cases resulted in severe population depression (%.1f%% avg utility loss).\n\n', ...
    overall_stats.Cases, overall_stats.UtilityLoss_pct);

fprintf('========================================\n');