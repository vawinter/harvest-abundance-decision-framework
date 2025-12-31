%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TRUE MISMATCH TESTING: Comparing Perceived vs Actual Population Trend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tests what happens when managers THINK the population is stable
% (use Stable weights + Stable parameters for decision-making)
% but the population is ACTUALLY decreasing (apply decision to Decrease parameters)
%
% Comparison:
% Run 1: Stable weights + Stable parameters → Optimal decision
% Run 2: Stable weights + Decrease parameters → See what happens with that decision
% Run 3: Decrease weights + Decrease parameters → What SHOULD be done
%
% Created by VAW
% Date: December 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

% Add MDPSolve path
addpath(genpath('C:/Users/vaw5154/OneDrive - The Pennsylvania State University/PhD/PSUTurkey/turkey_SDP/'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_results = true;

% Create output directory
output_dir = 'Results/TrueMismatch_Results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = (0:3)';
slope = 0.2;
uw_fixed = 0.2;
omu = 12.998;

% Simulation parameters
reps = 100000;
time_steps = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE COMBINATIONS TO TEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
months = {'January', 'April', 'September'};
scenarios = {'A', 'B', 'C'};
weathers = {'warm', 'cold'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STORAGE FOR RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_results
    all_results = struct();
    result_counter = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m = 1:length(months)
    scenario_month = months{m};
    
    % Weather loop (only for April)
    if strcmp(scenario_month, 'April')
        weather_list = weathers;
    else
        weather_list = {''};
    end
    
    for w = 1:length(weather_list)
        april_weather = weather_list{w};
        
        for s = 1:length(scenarios)
            scenario_weight = scenarios{s};
            
            % Create unique identifier
            if strcmp(scenario_month, 'April')
                run_id = sprintf('%s_%s_Scenario%s', scenario_month, april_weather, scenario_weight);
            else
                run_id = sprintf('%s_Scenario%s', scenario_month, scenario_weight);
            end
            
            fprintf('\n========================================\n');
            fprintf('Running: %s\n', run_id);
            fprintf('========================================\n');
            
            try
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% RUN 1: Manager's Decision (thinks population is Stable)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf('RUN 1: Manager thinks STABLE (making decision)...\n');
                
                % Set STABLE parameters
                [Fbar_stable, Pbar_stable, osig_stable, use_w_stable] = ...
                    setParameters(scenario_month, 'Stable', april_weather);
                
                % Load STABLE weights
                fname = sprintf('Data/norm_weights_popsize_scenario_%s.csv', scenario_weight);
                T = readtable(fname, 'VariableNamingRule', 'preserve');
                row = T(strcmp(T.("Population.Size"), 'Stable'), :);
                setGlobalWeights(row.normalized_weights);
                
                % Run model with STABLE parameters
                [model_stable, results_stable, mm_stable, ES_stable, ~, pp_stable, F_stable, ~, ~, ~, D_stable, ~, ~] = ...
                    PennTurkeyModel(uw_fixed, Fbar_stable, Pbar_stable, slope, osig_stable, omu, use_w_stable);
                
                % Get optimal decision from STABLE model
                if ~isempty(ES_stable) && length(ES_stable) >= 4
                    distances = sqrt((results_stable.Xopt(:,2) - ES_stable(2)).^2 + ...
                                   (results_stable.Xopt(:,3) - ES_stable(3)).^2 + ...
                                   (results_stable.Xopt(:,4) - ES_stable(4)).^2);
                    [~, closest_idx] = min(distances);
                    optimal_decision_stable = results_stable.Xopt(closest_idx, 1);
                else
                    optimal_decision_stable = mode(model_stable.X(results_stable.Ixopt, 1));
                end
                
                fprintf('  Manager chooses: %d weeks (based on Stable)\n', optimal_decision_stable);
                
                % Calculate utility under STABLE assumptions
                if ~isempty(ES_stable) && length(ES_stable) >= 4 && ~isempty(results_stable.v)
                    distances = sqrt((model_stable.X(:,2) - ES_stable(2)).^2 + ...
                                   (model_stable.X(:,3) - ES_stable(3)).^2 + ...
                                   (model_stable.X(:,4) - ES_stable(4)).^2);
                    [~, closest_idx] = min(distances);
                    utility_stable = results_stable.v(closest_idx);
                elseif ~isempty(pp_stable) && ~isempty(results_stable.v)
                    utility_stable = pp_stable' * results_stable.v;
                else
                    utility_stable = mean(results_stable.v);
                end
                
                % Simulate outcomes under STABLE
                [F_means_stable, MJ_means_stable] = runSimulations(D_stable, L, reps, time_steps);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% RUN 2: Reality Check (population is actually Decreasing)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf('RUN 2: Reality is DECREASE (applying decision)...\n');
                
                % Set DECREASE parameters
                [Fbar_decrease, Pbar_decrease, osig_decrease, use_w_decrease] = ...
                    setParameters(scenario_month, 'Decrease', april_weather);
                
                % Keep STABLE weights (manager still thinks it's stable)
                setGlobalWeights(row.normalized_weights);
                
                % Run model with DECREASE parameters
                [model_decrease, results_decrease, mm_decrease, ES_decrease, ~, pp_decrease, F_decrease, ~, ~, ~, D_decrease, ~, ~] = ...
                    PennTurkeyModel(uw_fixed, Fbar_decrease, Pbar_decrease, slope, osig_decrease, omu, use_w_decrease);
                
                % What WOULD be optimal under DECREASE?
                if ~isempty(ES_decrease) && length(ES_decrease) >= 4
                    distances = sqrt((results_decrease.Xopt(:,2) - ES_decrease(2)).^2 + ...
                                   (results_decrease.Xopt(:,3) - ES_decrease(3)).^2 + ...
                                   (results_decrease.Xopt(:,4) - ES_decrease(4)).^2);
                    [~, closest_idx] = min(distances);
                    optimal_decision_decrease = results_decrease.Xopt(closest_idx, 1);
                else
                    optimal_decision_decrease = mode(model_decrease.X(results_decrease.Ixopt, 1));
                end
                
                fprintf('  Reality optimal: %d weeks (if they knew)\n', optimal_decision_decrease);
                fprintf('  Decision error: %d weeks\n', optimal_decision_stable - optimal_decision_decrease);
                
                % Calculate utility under DECREASE
                if ~isempty(ES_decrease) && length(ES_decrease) >= 4 && ~isempty(results_decrease.v)
                    distances = sqrt((model_decrease.X(:,2) - ES_decrease(2)).^2 + ...
                                   (model_decrease.X(:,3) - ES_decrease(3)).^2 + ...
                                   (model_decrease.X(:,4) - ES_decrease(4)).^2);
                    [~, closest_idx] = min(distances);
                    utility_decrease = results_decrease.v(closest_idx);
                elseif ~isempty(pp_decrease) && ~isempty(results_decrease.v)
                    utility_decrease = pp_decrease' * results_decrease.v;
                else
                    utility_decrease = mean(results_decrease.v);
                end
                
                % Simulate outcomes under DECREASE
                [F_means_decrease, MJ_means_decrease] = runSimulations(D_decrease, L, reps, time_steps);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% RUN 3: Correct Decision (knows it's Decreasing)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf('RUN 3: Manager knows DECREASE (correct decision)...\n');
                
                % Load DECREASE weights (what they should use)
                row_decrease = T(strcmp(T.("Population.Size"), 'Decrease'), :);
                setGlobalWeights(row_decrease.normalized_weights);
                
                % Run model with DECREASE parameters and DECREASE weights
                [model_correct, results_correct, mm_correct, ES_correct, ~, pp_correct, ~, ~, ~, ~, D_correct, ~, ~] = ...
                    PennTurkeyModel(uw_fixed, Fbar_decrease, Pbar_decrease, slope, osig_decrease, omu, use_w_decrease);
                
                % Get optimal decision with correct knowledge
                if ~isempty(ES_correct) && length(ES_correct) >= 4
                    distances = sqrt((results_correct.Xopt(:,2) - ES_correct(2)).^2 + ...
                                   (results_correct.Xopt(:,3) - ES_correct(3)).^2 + ...
                                   (results_correct.Xopt(:,4) - ES_correct(4)).^2);
                    [~, closest_idx] = min(distances);
                    optimal_decision_correct = results_correct.Xopt(closest_idx, 1);
                else
                    optimal_decision_correct = mode(model_correct.X(results_correct.Ixopt, 1));
                end
                
                fprintf('  Correct decision: %d weeks (with Decrease weights)\n', optimal_decision_correct);
                
                % Calculate utility with correct knowledge
                if ~isempty(ES_correct) && length(ES_correct) >= 4 && ~isempty(results_correct.v)
                    distances = sqrt((model_correct.X(:,2) - ES_correct(2)).^2 + ...
                                   (model_correct.X(:,3) - ES_correct(3)).^2 + ...
                                   (model_correct.X(:,4) - ES_correct(4)).^2);
                    [~, closest_idx] = min(distances);
                    utility_correct = results_correct.v(closest_idx);
                elseif ~isempty(pp_correct) && ~isempty(results_correct.v)
                    utility_correct = pp_correct' * results_correct.v;
                else
                    utility_correct = mean(results_correct.v);
                end
                
                % Simulate outcomes with correct decision
                [F_means_correct, MJ_means_correct] = runSimulations(D_correct, L, reps, time_steps);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Calculate consequences of mismatch
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Population outcomes if manager's decision is applied to decreasing pop
                female_at_manager_decision = F_means_decrease(optimal_decision_stable + 1);
                male_at_manager_decision = MJ_means_decrease(optimal_decision_stable + 1);
                
                % Population outcomes if correct decision is made
                female_at_correct_decision = F_means_correct(optimal_decision_correct + 1);
                male_at_correct_decision = MJ_means_correct(optimal_decision_correct + 1);
                
                % Utility loss from mismatch
                utility_loss = utility_correct - utility_decrease;
                
                fprintf('\n  CONSEQUENCES:\n');
                fprintf('  Female: %.3f (manager) vs %.3f (correct) → Loss: %.3f\n', ...
                    female_at_manager_decision, female_at_correct_decision, ...
                    female_at_correct_decision - female_at_manager_decision);
                fprintf('  Utility loss from mismatch: %.3f\n', utility_loss);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Store results
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if save_results
                    all_results(result_counter).run_id = run_id;
                    all_results(result_counter).month = scenario_month;
                    all_results(result_counter).weather = april_weather;
                    all_results(result_counter).scenario = scenario_weight;
                    
                    % Manager's decision (thinks Stable)
                    all_results(result_counter).manager_decision = optimal_decision_stable;
                    all_results(result_counter).manager_utility_assumed = utility_stable;
                    all_results(result_counter).manager_female_assumed = F_means_stable(optimal_decision_stable + 1);
                    all_results(result_counter).manager_male_assumed = MJ_means_stable(optimal_decision_stable + 1);
                    
                    % Reality (actually Decrease, using manager's decision)
                    all_results(result_counter).reality_decision = optimal_decision_stable;
                    all_results(result_counter).reality_utility = utility_decrease;
                    all_results(result_counter).reality_female = female_at_manager_decision;
                    all_results(result_counter).reality_male = male_at_manager_decision;
                    
                    % Correct (knows Decrease)
                    all_results(result_counter).correct_decision = optimal_decision_correct;
                    all_results(result_counter).correct_utility = utility_correct;
                    all_results(result_counter).correct_female = female_at_correct_decision;
                    all_results(result_counter).correct_male = male_at_correct_decision;
                    
                    % Mismatch consequences
                    all_results(result_counter).decision_error = optimal_decision_stable - optimal_decision_correct;
                    all_results(result_counter).female_loss = female_at_correct_decision - female_at_manager_decision;
                    all_results(result_counter).male_loss = male_at_correct_decision - male_at_manager_decision;
                    all_results(result_counter).utility_loss = utility_loss;
                    
                    % Store all population outcomes for plotting
                    all_results(result_counter).F_stable = F_means_stable;
                    all_results(result_counter).MJ_stable = MJ_means_stable;
                    all_results(result_counter).F_decrease = F_means_decrease;
                    all_results(result_counter).MJ_decrease = MJ_means_decrease;
                    all_results(result_counter).F_correct = F_means_correct;
                    all_results(result_counter).MJ_correct = MJ_means_correct;
                    
                    result_counter = result_counter + 1;
                end
                
            catch ME
                warning('Error in run %s: %s', run_id, ME.message);
                fprintf('Error details: %s\n', ME.getReport());
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_results && exist('all_results', 'var') && ~isempty(all_results)
    % Save full structure
    save(fullfile(output_dir, 'true_mismatch_results.mat'), 'all_results');
    
    % Create CSV summary
    results_simple = struct();
    for i = 1:length(all_results)
        results_simple(i).run_id = all_results(i).run_id;
        results_simple(i).month = all_results(i).month;
        results_simple(i).weather = all_results(i).weather;
        results_simple(i).scenario = all_results(i).scenario;
        results_simple(i).manager_decision = all_results(i).manager_decision;
        results_simple(i).correct_decision = all_results(i).correct_decision;
        results_simple(i).decision_error = all_results(i).decision_error;
        results_simple(i).manager_female_assumed = all_results(i).manager_female_assumed;
        results_simple(i).reality_female = all_results(i).reality_female;
        results_simple(i).correct_female = all_results(i).correct_female;
        results_simple(i).female_loss = all_results(i).female_loss;
        results_simple(i).manager_utility_assumed = all_results(i).manager_utility_assumed;
        results_simple(i).reality_utility = all_results(i).reality_utility;
        results_simple(i).correct_utility = all_results(i).correct_utility;
        results_simple(i).utility_loss = all_results(i).utility_loss;
    end
    
    results_table = struct2table(results_simple);
    writetable(results_table, fullfile(output_dir, 'true_mismatch_summary.csv'));
    
    fprintf('\n========================================\n');
    fprintf('TRUE MISMATCH ANALYSIS COMPLETE\n');
    fprintf('Results saved: %d runs\n', length(all_results));
    fprintf('Files saved to: %s\n', output_dir);
    fprintf('========================================\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Fbar, Pbar, osig, use_w] = setParameters(scenario_month, scenario_trend, april_weather)
    if strcmp(scenario_trend, 'Increase')
        Fbar = 2.75;
        Pbar_base = 1.98;
    elseif strcmp(scenario_trend, 'Stable')
        Fbar = 2.5;
        Pbar_base = 1.8;
    else  % Decrease
        Fbar = 2.25;
        Pbar_base = 1.62;
    end
    
    if strcmp(scenario_month, 'April')
        if strcmp(april_weather, 'cold')
            Pbar = Pbar_base - 0.035;
        else
            Pbar = Pbar_base + 0.039;
        end
    else
        Pbar = Pbar_base;
    end
    
    if strcmp(scenario_month, 'September')
        osig = 1.0;
        use_w = false;
    else
        osig = 2.3;
        use_w = true;
    end
end

function setGlobalWeights(weights)
    global w0 w1 w2 w3
    w0 = weights(1);
    w1 = weights(2);
    w2 = weights(3);
    w3 = weights(4);
end

function [F_means, MJ_means] = runSimulations(D, L, reps, time_steps)
    F_means = zeros(length(L), 1);
    MJ_means = zeros(length(L), 1);
    
    for i = 1:length(L)
        SS = dsim(D, ones(reps,1) * [1.5, 1.5, 3], time_steps, L(i), [], [], 0);
        F_means(i) = mean(SS{4}(:, end));
        MJ_means(i) = mean(SS{1}(:, end) + SS{2}(:, end));
    end
end

function uw = getUtilityWeight(L)
    global w0 w1 w2 w3
    switch L
        case 1, uw = w1;
        case 2, uw = w2;
        case 3, uw = w3;
        otherwise, uw = w0;
    end
end

%% PennTurkeyModel function (same as before)
function [model, results, mm, ES, aa, pp, F, J, M, svals, D, osig, omu] = PennTurkeyModel(uw_fixed, Fbar, Pbar, slope, osig, omu, use_w)
    addfactor = 0.1;
    
    addfactor1_table = zeros(4,1);
    addfactor2_table = zeros(4,1);
    for i = 0:3
        hunter_pref = getUtilityWeight(i);
        if hunter_pref == 1
            addfactor1_table(i+1) = 0.1;
        else
            addfactor1_table(i+1) = 1/(2^(1/hunter_pref)-1);
        end
        if hunter_pref == 0
            addfactor2_table(i+1) = 0.1;
        else
            addfactor2_table(i+1) = 1/(2^(1/(1-hunter_pref))-1);
        end
    end
    addfactor2_table = zeros(4,1);
    
    svals = [];
    pp = [];
    aa = [];
    
    delta = 0.99;
    gammam = 0.409;
    gammaj = 0.653;
    gammaf = 0.542;
    cvv = 0.0736;
    cvw = 1.116 / 4.96;
    eta = 1 + 2 * Fbar * slope / Pbar;
    
    L = (0:3)';
    O = (0:26)';
    po = normpdf(O, omu, osig);
    po(O <= 5) = sum(po(O <= 5));
    po(O >= 15) = sum(po(O >= 15));
    po = po(6:16) / sum(po(6:16));
    H = [zeros(11,1) max(0.01,(0.07 - 0.00686*O(6:16) + 0.0175*(1:4)))];
    eH = H' * po;
    
    Ms = @(Md, v) gammam * Md .* v;
    Js = @(Jd, v) gammaj * Jd .* v;
    Fs = @(Fd, v) gammaf * Fd .* v;
    H = @(L) eH(L+1);
    
    if use_w
        Ps = @(F, w) 2 * eta * Pbar ./ (eta + 1 + (eta - 1) * (F / Fbar).^eta) .* w;
    else
        Ps = @(F, w) 2 * eta * Pbar ./ (eta + 1 + (eta - 1) * (F / Fbar).^eta) * 1;
    end
    
    Md = @(Ms, Js) Ms + Js;
    Jd = @(Fs, Pa, H) (1 - H) .* Fs .* Pa / 2;
    Fd = @(Fs, Pa, H) (1 - H) .* Fs .* (1 + Pa / 2);
    
    vparams = Burr3mom2param(1, cvv, 0);
    wparams = Burr3mom2param(1, cvw, 0);
    v = rvdef('burr3', vparams, 25);
    w = rvdef('burr3', wparams, 25);
    
    minpop = 0.00;
    M = linspace(minpop, 4, 41)';
    J = linspace(minpop, 4, 41)';
    F = linspace(minpop, 7, 71)';
    
    ut = @(L, Md, Jd, Fd) arrayfun(@(l, md, jd, fd) ...
        ifthenelse(fd == F(1) && l > 0, -inf, ...
        ((l + addfactor1_table(l+1)).^uw_fixed .* ...
         (md + jd + addfactor2_table(l+1)).^(1-uw_fixed)) .* (fd > F(1)) - l * 1e-10), ...
        L, Md, Jd, Fd);
    
    D = [];
    D = add2diagram(D, 'L', 'a', true, {}, L, [0.1278, 0.3713], []);
    D = add2diagram(D, 'Mj', 's', true, {}, M, [0.1196, 0.7941], []);
    D = add2diagram(D, 'Jj', 's', true, {}, J, [0.1125, 0.7143], []);
    D = add2diagram(D, 'Fj', 's', true, {}, F, [0.1172, 0.6160], []);
    D = add2diagram(D, 'v', 'c', true, {}, v, [0.3005, 0.9365], []);
    D = add2diagram(D, 'w', 'c', true, {}, w, [0.3282, 0.4451], []);
    D = add2diagram(D, 'H', 'c', true, {'L'}, H, [0.5181, 0.3684], [5 1]);
    D = add2diagram(D, 'Ms', 'c', true, {'Mj','v'}, Ms, [0.4804, 0.7865], [5 1; 5 1]);
    D = add2diagram(D, 'Js', 'c', true, {'Jj','v'}, Js, [0.4810, 0.7025], [5 1; 5 1]);
    D = add2diagram(D, 'Fs', 'c', true, {'Fj','v'}, Fs, [0.4751, 0.6156], [5 1; 5 1]);
    D = add2diagram(D, 'Ps', 'c', true, {'Fj','w'}, Ps, [0.4760, 0.4859], [5 1; 5 1]);
    D = add2diagram(D, 'Mj+', 'f', true, {'Ms','Js'}, Md, [0.7376, 0.7955], [5 1; 5 1]);
    D = add2diagram(D, 'Jj+', 'f', true, {'Fs','Ps','H'}, Jd, [0.7477, 0.7058], [5 1; 5 1; 5 1]);
    D = add2diagram(D, 'Fj+', 'f', true, {'Fs','Ps','H'}, Fd, [0.7706, 0.6217], [5 1; 5 1; 5 1]);
    D = add2diagram(D, 'utility', 'u', true, {'L','Mj','Jj','Fj'}, ut, [0.3649, 0.1564], [5 1; 5 1; 5 1; 5 1]);
    
    doptions = struct('d', delta, 'cleanup', 2, 'reps', 0, 'chunk', 1000, 'print', 1, 'ptype', 0);
    model = d2model(D, doptions);
    
    moptions = struct('print', 0);
    results = mdpsolve(model, moptions);
    
    if isfield(results, 'Ixopt')
        Xopt = model.X(results.Ixopt, :);
        results.Xopt = Xopt;
    else
        results.Xopt = [];
    end
    
    if ~isempty(results.pstar)
        ns=size(results.pstar,1); pp=ones(ns,1)/ns;
        for i=1:500, pp=results.pstar*pp; end
        mm=marginals(pp,[length(M) length(J) length(F)]);
        ES=pp'*[results.Xopt Ps(results.Xopt(:,4),1)];
        if nargout>=3
            aa=zeros(1,4);
            for i=0:3
                aa(i+1)=sum(pp(results.Xopt(:,1)==i));
            end
        end
    else
        reps=100000;
        SS=dsim(D,ones(reps,1)*[1.5 1.5 3],200,model.X(results.Ixopt,1),[],[],0);
        svals=linspace(0,max(F),501)';
        mm=cell(1,5);
        mm{1}=histcounts(SS{1},[(0:4)';inf])'/reps;
        for i=2:4
            mm{i}=kernel(svals,SS{i});
        end
        mm{5}=kernel(svals,SS{2}+SS{3});
        ES = [mean(SS{1}) mean(SS{2}) mean(SS{3}) mean(SS{2}+SS{3}) mean(SS{4})];
    end
end