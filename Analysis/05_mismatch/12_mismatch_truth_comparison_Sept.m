%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TRUE MISMATCH TESTING: Comparing Perceived vs Actual Population Trend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests what happens when managers THINK the population is Stable
% but the population is ACTUALLY Decreasing.
%
% Run 1: Stable weights + Stable parameters     → Manager's assumed optimal
% Run 2: Stable weights + Decrease parameters   → Reality under manager's decision
% Run 3: Decrease weights + Decrease parameters → What SHOULD be done
%
% September decision point: mast observed, recruitment known (cvw=0, osig=1.0)
%
% Winter et al. 20XX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

addpath(genpath('C:/Users/vaw5154/OneDrive - The Pennsylvania State University/PhD/PSUTurkey/turkey_SDP/'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
slope    = 0.2;
uw_fixed = 0.2;
omu      = 12.998;   % change to 5 or 16 for mast sensitivity
L        = (0:3)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMBINATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
months    = {'September'};
scenarios = {'A', 'B', 'C'};
weathers  = {'warm', 'cold'};

if save_results
    all_results    = struct();
    result_counter = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:length(months)
    scenario_month = months{m};

    % September and April both use weather to set Pbar
    % January does not
    if strcmp(scenario_month, 'January')
        weather_list = {''};
    else
        weather_list = weathers;
    end

    for w = 1:length(weather_list)
        april_weather = weather_list{w};

        for s = 1:length(scenarios)
            scenario_weight = scenarios{s};

            % Run identifier
            if strcmp(scenario_month, 'January')
                run_id = sprintf('%s_Scenario%s', scenario_month, scenario_weight);
            else
                run_id = sprintf('%s_%s_Scenario%s', scenario_month, april_weather, scenario_weight);
            end

            fprintf('\n========================================\n');
            fprintf('Running: %s\n', run_id);
            fprintf('========================================\n');

            try
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% LOAD WEIGHT TABLE ONCE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fname = sprintf('Data/norm_weights_popsize_scenario_%s.csv', scenario_weight);
                T     = readtable(fname, 'VariableNamingRule', 'preserve');

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% RUN 1: Manager thinks Stable -- makes decision
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf('RUN 1: Manager thinks STABLE...\n');

                [Fbar_s, Pbar_s, osig_s, cvw_s] = setParameters(scenario_month, 'Stable', april_weather);
                use_w_s = (cvw_s > 0);

                row_stable = T(strcmp(T.("Population.Size"), 'Stable'), :);
                setGlobalWeights(row_stable.normalized_weights);

                [~, results_stable, ~, ES_stable, ~, pp_stable, ~, ~, ~, ~, ~, ~, ~] = ...
                    PennTurkeyModel(uw_fixed, Fbar_s, Pbar_s, slope, osig_s, omu, use_w_s, cvw_s);

                optimal_stable   = getOptimalDecision(results_stable, ES_stable);
                utility_stable   = pp_stable' * results_stable.v;
                [F_means_s, MJ_means_s] = computePopResponse(L, results_stable, pp_stable);

                fprintf('  Manager chooses: %d weeks\n', optimal_stable);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% RUN 2: Reality is Decreasing -- apply manager's decision
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf('RUN 2: Reality is DECREASE (manager unaware)...\n');

                [Fbar_d, Pbar_d, osig_d, cvw_d] = setParameters(scenario_month, 'Decrease', april_weather);
                use_w_d = (cvw_d > 0);

                % Manager still uses Stable weights -- unaware of true trend
                setGlobalWeights(row_stable.normalized_weights);

                [~, results_decrease, ~, ES_decrease, ~, pp_decrease, ~, ~, ~, ~, ~, ~, ~] = ...
                    PennTurkeyModel(uw_fixed, Fbar_d, Pbar_d, slope, osig_d, omu, use_w_d, cvw_d);

                optimal_decrease  = getOptimalDecision(results_decrease, ES_decrease);
                utility_decrease  = pp_decrease' * results_decrease.v;
                [F_means_d, MJ_means_d] = computePopResponse(L, results_decrease, pp_decrease);

                fprintf('  Reality optimal (if known): %d weeks\n', optimal_decrease);
                fprintf('  Decision error: %d weeks\n', optimal_stable - optimal_decrease);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% RUN 3: Correct -- manager knows it is Decreasing
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf('RUN 3: Manager knows DECREASE -- correct decision...\n');

                row_decrease = T(strcmp(T.("Population.Size"), 'Decrease'), :);
                setGlobalWeights(row_decrease.normalized_weights);

                [~, results_correct, ~, ES_correct, ~, pp_correct, ~, ~, ~, ~, ~, ~, ~] = ...
                    PennTurkeyModel(uw_fixed, Fbar_d, Pbar_d, slope, osig_d, omu, use_w_d, cvw_d);

                optimal_correct  = getOptimalDecision(results_correct, ES_correct);
                utility_correct  = pp_correct' * results_correct.v;
                [F_means_c, MJ_means_c] = computePopResponse(L, results_correct, pp_correct);

                fprintf('  Correct decision: %d weeks\n', optimal_correct);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% CONSEQUENCES OF MISMATCH
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Population outcomes when manager's decision (optimal_stable)
                % is applied to the decreasing population
                female_manager  = F_means_d(optimal_stable + 1);
                male_manager    = MJ_means_d(optimal_stable + 1);

                % Population outcomes under correct decision
                female_correct  = F_means_c(optimal_correct + 1);
                male_correct    = MJ_means_c(optimal_correct + 1);

                % Utility loss = what you would have gotten - what you got
                utility_loss = utility_correct - utility_decrease;

                fprintf('\n  CONSEQUENCES:\n');
                fprintf('  Female: %.3f (manager) vs %.3f (correct) -- Loss: %.3f\n', ...
                    female_manager, female_correct, female_correct - female_manager);
                fprintf('  Utility loss from mismatch: %.4f\n', utility_loss);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% STORE RESULTS
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if save_results
                    r = result_counter;

                    all_results(r).run_id          = run_id;
                    all_results(r).month           = scenario_month;
                    all_results(r).weather         = april_weather;
                    all_results(r).scenario        = scenario_weight;

                    % Manager (thinks Stable)
                    all_results(r).manager_decision        = optimal_stable;
                    all_results(r).manager_utility_assumed = utility_stable;
                    all_results(r).manager_female_assumed  = F_means_s(optimal_stable + 1);
                    all_results(r).manager_male_assumed    = MJ_means_s(optimal_stable + 1);

                    % Reality (Decrease, using manager's decision)
                    all_results(r).reality_decision = optimal_stable;
                    all_results(r).reality_utility  = utility_decrease;
                    all_results(r).reality_female   = female_manager;
                    all_results(r).reality_male     = male_manager;

                    % Correct (knows Decrease)
                    all_results(r).correct_decision = optimal_correct;
                    all_results(r).correct_utility  = utility_correct;
                    all_results(r).correct_female   = female_correct;
                    all_results(r).correct_male     = male_correct;

                    % Mismatch consequences
                    all_results(r).decision_error = optimal_stable - optimal_correct;
                    all_results(r).female_loss    = female_correct - female_manager;
                    all_results(r).male_loss      = male_correct   - male_manager;
                    all_results(r).utility_loss   = utility_loss;

                    % Full population response curves for plotting
                    all_results(r).F_stable   = F_means_s;
                    all_results(r).MJ_stable  = MJ_means_s;
                    all_results(r).F_decrease = F_means_d;
                    all_results(r).MJ_decrease = MJ_means_d;
                    all_results(r).F_correct  = F_means_c;
                    all_results(r).MJ_correct = MJ_means_c;

                    result_counter = result_counter + 1;
                end

            catch ME
                warning('Error in run %s: %s', run_id, ME.message);
                fprintf('Error details:\n%s\n', ME.getReport());
            end

        end  % scenarios
    end  % weather
end  % months

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_results && exist('all_results', 'var') && ~isempty(all_results)
    save(fullfile(output_dir, 'true_mismatch_results.mat'), 'all_results');

    % Scalar summary for CSV
    results_simple = struct();
    for i = 1:length(all_results)
        results_simple(i).run_id                 = all_results(i).run_id;
        results_simple(i).month                  = all_results(i).month;
        results_simple(i).weather                = all_results(i).weather;
        results_simple(i).scenario               = all_results(i).scenario;
        results_simple(i).manager_decision       = all_results(i).manager_decision;
        results_simple(i).correct_decision       = all_results(i).correct_decision;
        results_simple(i).decision_error         = all_results(i).decision_error;
        results_simple(i).manager_female_assumed = all_results(i).manager_female_assumed;
        results_simple(i).reality_female         = all_results(i).reality_female;
        results_simple(i).correct_female         = all_results(i).correct_female;
        results_simple(i).female_loss            = all_results(i).female_loss;
        results_simple(i).manager_utility_assumed = all_results(i).manager_utility_assumed;
        results_simple(i).reality_utility        = all_results(i).reality_utility;
        results_simple(i).correct_utility        = all_results(i).correct_utility;
        results_simple(i).utility_loss           = all_results(i).utility_loss;
    end

    writetable(struct2table(results_simple), fullfile(output_dir, 'true_mismatch_summary.csv'));

    fprintf('\n========================================\n');
    fprintf('MISMATCH ANALYSIS COMPLETE\n');
    fprintf('Results saved: %d runs to %s\n', length(all_results), output_dir);
    fprintf('========================================\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Fbar, Pbar, osig, cvw] = setParameters(scenario_month, scenario_trend, april_weather)
    switch scenario_trend
        case 'Increase', Fbar = 2.75; Pbar_base = 1.98;
        case 'Stable',   Fbar = 2.50; Pbar_base = 1.80;
        otherwise,       Fbar = 2.00; Pbar_base = 1.40;  % Decrease
    end

    if strcmp(scenario_month, 'April') || strcmp(scenario_month, 'September')
        switch april_weather
            case 'cold', Pbar = Pbar_base - 0.035;
            case 'warm', Pbar = Pbar_base + 0.039;
            otherwise,   Pbar = Pbar_base;
        end
    else
        Pbar = Pbar_base;
    end

    switch scenario_month
        case 'September'
            osig = 1.0;
            cvw  = 0;
        case 'April'
            osig = 2.3;
            cvw  = 1.116/4.96;
        otherwise  % January
            osig = 2.3;
            cvw  = 1.116/4.96;
    end
end

%--------------------------------------------------------------------------
function setGlobalWeights(weights)
    global w0 w1 w2 w3
    w0 = weights(1);
    w1 = weights(2);
    w2 = weights(3);
    w3 = weights(4);
end

%--------------------------------------------------------------------------
function uw = getUtilityWeight(L)
    global w0 w1 w2 w3
    switch L
        case 1,    uw = w1;
        case 2,    uw = w2;
        case 3,    uw = w3;
        otherwise, uw = w0;
    end
end

%--------------------------------------------------------------------------
function opt = getOptimalDecision(results, ES)
    % Returns optimal season length at the long-run expected state
    if ~isempty(ES) && length(ES) >= 4
        distances = sqrt((results.Xopt(:,2) - ES(2)).^2 + ...
                         (results.Xopt(:,3) - ES(3)).^2 + ...
                         (results.Xopt(:,4) - ES(4)).^2);
        [~, idx] = min(distances);
        opt = results.Xopt(idx, 1);
    else
        opt = mode(results.Xopt(:, 1));
    end
end

%--------------------------------------------------------------------------
function [F_means, MJ_means] = computePopResponse(L, results, pp)
    % Stationary-distribution-weighted mean populations by season length
    F_means  = zeros(length(L), 1);
    MJ_means = zeros(length(L), 1);
    for i = 1:length(L)
        idx = results.Xopt(:,1) == L(i);
        w   = pp(idx);
        sw  = sum(w);
        if sw > 0
            F_means(i)  = w' * results.Xopt(idx, 4) / sw;
            MJ_means(i) = w' * (results.Xopt(idx,2) + results.Xopt(idx,3)) / sw;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DYNAMIC MODEL FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model, results, mm, ES, aa, pp, F, J, M, svals, D, osig, omu] = ...
        PennTurkeyModel(uw_fixed, Fbar, Pbar, slope, osig, omu, use_w, cvw)

    %% Utility weight addfactor tables (season lengths 0-3)
    addfactor1_table = zeros(4,1);
    addfactor2_table = zeros(4,1);
    for i = 0:3
        hp = getUtilityWeight(i);
        addfactor1_table(i+1) = ifthenelse(hp == 1, 0.1, 1/(2^(1/hp)-1));
        addfactor2_table(i+1) = 0;
    end

    svals = [];
    pp    = [];
    aa    = [];

    %% Population model parameters
    delta  = 0.99;
    gammam = 0.409;
    gammaj = 0.653;
    gammaf = 0.542;
    cvv    = 0.0736;

    eta = 1 + 2*Fbar*slope/Pbar;

    %% Decision variable
    L = (0:3)';

    %% Mast distribution (truncated normal, discretized over 5:15)
    Ovec = (0:26)';
    po   = pdfn(Ovec, omu, osig);
    po   = ifthenelse(Ovec<=5,  sum(po(Ovec<=5)),  ...
           ifthenelse(Ovec>=15, sum(po(Ovec>=15)), po));
    po   = po(6:16) / sum(po(6:16));
    Ogrid = Ovec(6:16);

    % Mast as i.i.d. discrete chance node (no parents, no transition)
    Odist.type       = 'discrete';
    Odist.parameters = [Ogrid, po];
    Odist.values     = Ogrid;
    Odist.cpt        = po;
    Odist.size       = length(Ogrid);
    Odist.lb         = min(Ogrid);
    Odist.simfunc    = [];
    Odist.ztype      = 0;

    %% Transition functions
    Ms = @(Md, v)          gammam * Md .* v;
    Js = @(Jd, v)          gammaj * Jd .* v;
    Fs = @(Fd, v)          gammaf * Fd .* v;
    H  = @(Lval, Oval)     max(0.01, 0.07 - 0.00686.*Oval + 0.0175.*Lval);

    if use_w
        Ps = @(F_val, w) 2*eta*Pbar ./ (eta+1+(eta-1)*(F_val/Fbar).^eta) .* w;
    else
        Ps = @(F_val, w) 2*eta*Pbar ./ (eta+1+(eta-1)*(F_val/Fbar).^eta);
    end

    Md = @(Ms_val, Js_val)        Ms_val + Js_val;
    Jd = @(Fs_val, Pa, H_val)     (1-H_val) .* Fs_val .* Pa / 2;
    Fd = @(Fs_val, Pa, H_val)     (1-H_val) .* Fs_val .* (1 + Pa/2);

    %% Random variables
    vparams = Burr3mom2param(1, cvv, 0);
    v = rvdef('burr3', vparams, 25);

    wparams = Burr3mom2param(1, max(cvw, 1e-6), 0);
    w = rvdef('burr3', wparams, 25);

    %% State variable grids
    minpop = 0.00;
    M = linspace(minpop, 4, 41)';
    J = linspace(minpop, 4, 41)';
    F = linspace(minpop, 7, 71)';

    %% Utility function
    ut = @(L_val, Md_val, Jd_val, Fd_val) arrayfun(@(l, md, jd, fd) ...
        ifthenelse(fd==F(1) && l>0, -inf, ...
        ((l + addfactor1_table(l+1)).^uw_fixed .* ...
         (md + jd + addfactor2_table(l+1)).^(1-uw_fixed)) .* (fd>F(1)) - l*1e-10), ...
        L_val, Md_val, Jd_val, Fd_val);

    %% Build decision diagram
    D = [];
    D = add2diagram(D, 'L',       'a', true, {},                  L,     [0.1278, 0.3713], []);
    D = add2diagram(D, 'Mj',      's', true, {},                  M,     [0.1196, 0.7941], []);
    D = add2diagram(D, 'Jj',      's', true, {},                  J,     [0.1125, 0.7143], []);
    D = add2diagram(D, 'Fj',      's', true, {},                  F,     [0.1172, 0.6160], []);
    D = add2diagram(D, 'O',       'c', true, {},                  Odist, [0.1100, 0.4500], []);
    D = add2diagram(D, 'v',       'c', true, {},                  v,     [0.3005, 0.9365], []);
    D = add2diagram(D, 'w',       'c', true, {},                  w,     [0.3282, 0.4451], []);
    D = add2diagram(D, 'H',       'c', true, {'L','O'},           H,     [0.5181, 0.3684], [length(L) 1; length(Ogrid) 1]);
    D = add2diagram(D, 'Ms',      'c', true, {'Mj','v'},          Ms,    [0.4804, 0.7865], [5 1; 5 1]);
    D = add2diagram(D, 'Js',      'c', true, {'Jj','v'},          Js,    [0.4810, 0.7025], [5 1; 5 1]);
    D = add2diagram(D, 'Fs',      'c', true, {'Fj','v'},          Fs,    [0.4751, 0.6156], [5 1; 5 1]);
    D = add2diagram(D, 'Ps',      'c', true, {'Fj','w'},          Ps,    [0.4760, 0.4859], [5 1; 5 1]);
    D = add2diagram(D, 'Mj+',     'f', true, {'Ms','Js'},         Md,    [0.7376, 0.7955], [5 1; 5 1]);
    D = add2diagram(D, 'Jj+',     'f', true, {'Fs','Ps','H'},     Jd,    [0.7477, 0.7058], [5 1; 5 1; 5 1]);
    D = add2diagram(D, 'Fj+',     'f', true, {'Fs','Ps','H'},     Fd,    [0.7706, 0.6217], [5 1; 5 1; 5 1]);
    D = add2diagram(D, 'utility', 'u', true, {'L','Mj','Jj','Fj'},ut,    [0.3649, 0.1564], [5 1; 5 1; 5 1; 5 1]);

    %% Solve MDP (ptype=1 ensures pstar is always computed)
    doptions = struct('d', delta, 'cleanup', 2, 'reps', 0, 'chunk', 1000, 'print', 1, 'ptype', 1);
    model    = d2model(D, doptions);

    moptions = struct('print', 0);
    results  = mdpsolve(model, moptions);

    fprintf('  pstar empty: %d, Ixopt size: %d\n', isempty(results.pstar), length(results.Ixopt));

    if isfield(results, 'Ixopt')
        results.Xopt = model.X(results.Ixopt, :);
    else
        results.Xopt = [];
    end

    %% Post-solution processing
    if ~isempty(results.pstar)
        ns = size(results.pstar, 1);
        pp = ones(ns,1) / ns;
        for i = 1:500
            pp = results.pstar * pp;
        end

        % Marginals over M, J, F (O integrated out by d2model)
        mm = marginals(pp, [length(M) length(J) length(F)]);
        ES = pp' * [results.Xopt, Ps(results.Xopt(:,4), 1)];

        % Season length distribution
        aa = zeros(4, 1);
        for i = 0:3
            aa(i+1) = sum(pp(results.Xopt(:,1) == i));
        end

    else
        %% Simulation fallback (should not be reached with ptype=1)
        warning('pstar empty -- falling back to simulation.');
        reps_sim = 100000;
        SS = dsim(D, ones(reps_sim,1) * [1.5, 1.5, 3, omu], 200, ...
                  model.X(results.Ixopt, 1), [], [], 0);
        svals = linspace(0, max(F), 501)';
        mm    = cell(1, 3);
        [mm{1}, ~] = ksdensity(SS{2}(:,end), svals);
        [mm{2}, ~] = ksdensity(SS{3}(:,end), svals);
        [mm{3}, ~] = ksdensity(SS{4}(:,end), svals);
        aa = zeros(4, 1);
        for i = 0:3
            aa(i+1) = mean(SS{1}(:,end) == i);
        end
        ES = [mean(SS{1}(:,end)) mean(SS{2}(:,end)) mean(SS{3}(:,end)) ...
              mean(SS{2}(:,end)+SS{3}(:,end)) mean(SS{4}(:,end))];
    end

end