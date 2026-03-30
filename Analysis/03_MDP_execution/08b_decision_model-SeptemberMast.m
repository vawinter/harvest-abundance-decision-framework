%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%X
%% A Decision Framework for Balancing Hunting Opportunity and Population Abundance: 
% A Case Study for Wild Turkey Management
% Winter et al. 20XX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%X
% Decision made at September time-point
% September = known recruitment and mast 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Automated Decision Framework Analysis - September
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

% Add MDPSolve path
addpath(genpath('C:/Users/vaw5154/OneDrive - The Pennsylvania State University/PhD/PSUTurkey/turkey_SDP/'))

% **USER OPTIONS**
save_plots = true;
save_results = true;

% Create output directory
output_dir = 'Results/Sept_Mast_StateVar_avg';
% output_dir = 'Results/Sept_Mast_StateVar_low';
% output_dir = 'Results/Sept_Mast_StateVar_high';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define all combinations to test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
months = {'September'};
trends = {
    %'Increase', 
'Stable'
%, 'Decrease'
};
scenarios = {'A', 'B', 'C'};
weathers = {'warm'
%, 'cold'
};
model_var = {'Low', 'High'};

% Fixed parameters
% Average mast
omu = 12.998;

% Low mast  
% omu = 6;

% High mast
% omu = 16;

slope  = 0.2;
uw_fixed = 0.2;
L = (0:3)';

% Simulation parameters
reps = 100000;
time_steps = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_results
    all_results = struct('run_id', {}, 'month', {}, 'trend', {}, 'scenario', {}, ...
                         'weather', {}, 'model_var', {}, 'value_function', {}, ...
                         'stationary_dist', {}, 'expected_utility', {}, ...
                         'optimal_action', {}, 'F_means', {}, 'MJ_means', {}, 'ES', {});
    result_counter = 1;
end

for m = 1:length(months)
    scenario_month = months{m};

    if strcmp(scenario_month, 'April')
        weather_list = weathers;
        var_list = model_var;
    elseif strcmp(scenario_month, 'September')
        weather_list = weathers;
        var_list = {''};
    else
        weather_list = {''};
        var_list = {''};
    end

    for v = 1:length(var_list)
        april_var = var_list{v};

        for w = 1:length(weather_list)
            april_weather = weather_list{w};

            for t = 1:length(trends)
                scenario_trend = trends{t};

                for s = 1:length(scenarios)
                    scenario_0 = scenarios{s};

                    if strcmp(scenario_month, 'April')
                        run_id = sprintf('%s_%s_%s_%s_%s', scenario_month, april_weather, april_var, scenario_trend, scenario_0);
                    elseif strcmp(scenario_month, 'September')
                        run_id = sprintf('%s_%s_%s_%s', scenario_month, april_weather, scenario_trend, scenario_0);
                    else
                        run_id = sprintf('%s_%s_%s', scenario_month, scenario_trend, scenario_0);
                    end

                    fprintf('\n========================================\n');
                    fprintf('Running: %s\n', run_id);
                    fprintf('========================================\n');

                    try
                        [Fbar, Pbar, osig, cvw] = setParameters(scenario_month, scenario_trend, april_weather, april_var);

                        fname = sprintf('Data/norm_weights_popsize_scenario_%s.csv', scenario_0);
                        T = readtable(fname, 'VariableNamingRule', 'preserve');
                        row = T(strcmp(T.("Population.Size"), scenario_trend), :);
                        setGlobalWeights(row.normalized_weights);

                        use_w = (cvw > 0);

                        [model, results, mm, ES, aa, pp, F, J, M, svals, D, ~, ~] = ...
                            PennTurkeyModel(uw_fixed, Fbar, Pbar, slope, osig, omu, use_w, cvw);

                        O = (6:16)';
                        O_mid = O(max(1, min(length(O), round((omu - O(1)) / (O(2) - O(1))) + 1)));
                        [F_means, MJ_means] = runSimulations(D, L, reps, time_steps, O_mid);

                       if save_results
                            all_results(result_counter).run_id    = run_id;
                            all_results(result_counter).month     = scenario_month;
                            all_results(result_counter).trend     = scenario_trend;
                            all_results(result_counter).scenario  = scenario_0;
                            all_results(result_counter).weather   = april_weather;
                            all_results(result_counter).value_function  = results.v;
                            all_results(result_counter).stationary_dist = pp;
                            all_results(result_counter).model_var = april_var;

                            
                            % Expected value function: E_π[V(s)] = pp' * v
                            if ~isempty(pp) && ~isempty(results.v)
                                % Stationary distribution available - exact EVF
                                all_results(result_counter).expected_utility = pp' * results.v;
                                all_results(result_counter).stationary_dist  = pp;
                            elseif ~isempty(results.v)
                                % Simulation fallback: weight v by empirical state visit frequencies
                                fprintf('Computing EVF via simulation weighting...\n');
                                reps_evf = 50000;
                                SS_evf = dsim(D, ones(reps_evf,1) * [1.5, 1.5, 3, O_mid], ...
                                            200, model.X(results.Ixopt, 1), [], [], 0);

                                M_step = M(2) - M(1);
                                J_step = J(2) - J(1);
                                F_step = F(2) - F(1);

                                M_idx = max(1, min(length(M), round((SS_evf{2}(:,end) - M(1)) / M_step) + 1));
                                J_idx = max(1, min(length(J), round((SS_evf{3}(:,end) - J(1)) / J_step) + 1));
                                F_idx = max(1, min(length(F), round((SS_evf{4}(:,end) - F(1)) / F_step) + 1));

                                nM = length(M); nJ = length(J); nF = length(F); nO = length(O);
                                v_grid = reshape(results.v, [nM nJ nF nO]);
                                v_MJF  = mean(v_grid, 4);

                                linear_idx = sub2ind([nM nJ nF], M_idx, J_idx, F_idx);

                                % Empirical state frequencies = approximation of π(s)
                                % Save as stationary_dist so VOI script can use for proper weighting
                                freq = histcounts(linear_idx, 1:nM*nJ*nF+1)' / reps_evf;

                                all_results(result_counter).expected_utility = freq' * v_MJF(:);
                                all_results(result_counter).stationary_dist  = freq;
                            else
                                warning('Could not calculate EVF for %s', run_id);
                                all_results(result_counter).expected_utility = NaN;
                                all_results(result_counter).stationary_dist  = [];
                            end

                            all_results(result_counter).optimal_action = mode(model.X(results.Ixopt, 1));
                            all_results(result_counter).F_means  = F_means;
                            all_results(result_counter).MJ_means = MJ_means;
                            all_results(result_counter).ES       = ES;
                            result_counter = result_counter + 1;
                        end

                        if save_plots
                            savePlots(results, aa, mm, ES, F, M, J, L, ...
                                      F_means, MJ_means, svals, run_id, output_dir);
                        end

                    catch ME
                        warning('Error in run %s: %s', run_id, ME.message);
                        fprintf('Error details: %s\n', ME.getReport());
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save resultsif 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_results && exist('all_results', 'var') && ...
        ~isempty(all_results) && isfield(all_results, 'run_id')

    save(fullfile(output_dir, 'all_utility_results.mat'), 'all_results');

    results_simple = struct();
    for i = 1:length(all_results)
        results_simple(i).run_id           = all_results(i).run_id;
        results_simple(i).month            = all_results(i).month;
        results_simple(i).trend            = all_results(i).trend;
        results_simple(i).scenario         = all_results(i).scenario;
        results_simple(i).weather          = all_results(i).weather;
        results_simple(i).model_var        = all_results(i).model_var;
        results_simple(i).expected_utility = all_results(i).expected_utility;
        results_simple(i).optimal_action   = all_results(i).optimal_action;
    end

    results_table = struct2table(results_simple);
    writetable(results_table, fullfile(output_dir, 'utility_results.csv'));

    fprintf('\n========================================\n');
    fprintf('Results saved: %d runs\n', length(all_results));
    fprintf('Files saved to: %s\n', output_dir);
    fprintf('========================================\n');
else
    fprintf('\n========================================\n');
    fprintf('No results to save\n');
    fprintf('========================================\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Fbar, Pbar, osig, cvw] = setParameters(scenario_month, scenario_trend, april_weather, april_var)
    if strcmp(scenario_trend, 'Increase')
        Fbar = 2.75;
        Pbar_base = 1.98;
    elseif strcmp(scenario_trend, 'Stable')
        Fbar = 2.5;
        Pbar_base = 1.8;
    else
        Fbar = 2.25;
        Pbar_base = 1.62;
    end

    if strcmp(scenario_month, 'April') || strcmp(scenario_month, 'September')
        if strcmp(april_weather, 'cold')
            Pbar = Pbar_base - 0.035;
        elseif strcmp(april_weather, 'warm')
            Pbar = Pbar_base + 0.039;
        else
            Pbar = Pbar_base;
        end
    else
        Pbar = Pbar_base;
    end

    if strcmp(scenario_month, 'September')
        osig = 1.0;
        cvw  = 0;
    elseif strcmp(scenario_month, 'April')
        osig = 2.3;
        if strcmp(april_var, 'Low')
            cvw = (1.116 / 4.96) * 0.5;
        else
            cvw = (1.116 / 4.96) * 2;
        end
    else
        osig = 2.3;
        cvw  = 1.116 / 4.96;
    end
end

function setGlobalWeights(weights)
    global w0 w1 w2 w3
    w0 = weights(1);
    w1 = weights(2);
    w2 = weights(3);
    w3 = weights(4);
end

function [F_means, MJ_means] = runSimulations(D, L, reps, time_steps, O_mid)
    F_means  = zeros(length(L), 1);
    MJ_means = zeros(length(L), 1);
    for i = 1:length(L)
        SS = dsim(D, ones(reps,1) * [1.5, 1.5, 3, O_mid], time_steps, L(i), [], [], 0);
        F_means(i)  = mean(SS{4}(:, end));
        MJ_means(i) = mean(SS{1}(:, end) + SS{2}(:, end));
    end
end

function savePlots(results, aa, mm, ES, F, M, J, L, F_means, MJ_means, ...
                   svals, run_id, output_dir)

    %% 1. Action plot
    if ~isempty(results.Xopt)
        fig = figure('Visible', 'off');
        try
            actionplot(results.Xopt, F, ES, 1);
            grid on;
            saveas(fig, fullfile(output_dir, sprintf('action_%s.png', run_id)));
        catch err
            warning('Could not create action plot for %s: %s', run_id, err.message);
        end
        close(fig);
    end

    %% 2. Long-run season length distribution
    if ~isempty(aa)
        fig = figure('Visible', 'off');
        L_plot = (0:length(aa)-1)';
        plot(L_plot, aa, 'k.-', 'LineWidth', 2, 'MarkerSize', 20);
        ylim([0 max(1, max(aa) * 1.1)]);
        xlim([-0.5, length(aa) - 0.5]);
        set(gca, 'XTick', 0:1:length(aa)-1, 'FontSize', 16);
        xlabel('Season Length (weeks)', 'Interpreter', 'latex', 'FontSize', 18);
        ylabel('Probability', 'Interpreter', 'latex', 'FontSize', 18);
      %  title(strrep(run_id, '_', '\_'), 'FontSize', 11);
        grid on;
        saveas(fig, fullfile(output_dir, sprintf('season_dist_%s.png', run_id)));
        close(fig);
    end

    %% 3. Long-run population density distributions
    if ~isempty(mm) && length(mm) >= 3
        fig = figure('Visible', 'off');
        set(fig, 'Position', [100, 100, 700, 500]);
        hold on;

        if length(mm{1}) == 501 && ~isempty(svals)
            xaxis_M = svals;
            xaxis_J = svals;
            xaxis_F = svals;
        else
            xaxis_M = M;
            xaxis_J = J;
            xaxis_F = F;
        end

        plot(xaxis_M, mm{1}, 'b-', 'LineWidth', 2, 'DisplayName', '$M$');
        plot(xaxis_J, mm{2}, 'r-', 'LineWidth', 2, 'DisplayName', '$J$');
        plot(xaxis_F, mm{3}, 'y-', 'LineWidth', 2, 'DisplayName', '$F$');

        hold off;
        legend('Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 14);
        xlim([0, max(F)]);
       % ylim([0, max([max(mm{1}), max(mm{2}), max(mm{3})]) * 1.15]);
       ylim([0, 3.5]);
        xlabel('Population Density', 'Interpreter', 'latex', 'FontSize', 16);
        ylabel('Probability Density', 'Interpreter', 'latex', 'FontSize', 16);
       % title(strrep(run_id, '_', '\_'), 'FontSize', 11);
        set(gca, 'FontSize', 14);
        grid on;
        saveas(fig, fullfile(output_dir, sprintf('longrun_pop_%s.png', run_id)));
        close(fig);
    end

    %% 4. Population response to season length
    fig = figure('Visible', 'off');
    set(fig, 'Position', [100, 100, 1200, 400]);
    subplot(1,2,1);
    plot(L, F_means, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Female Population by Season Length');
    xlabel('Season Length (weeks)');
    ylabel('Mean Female Density');
    grid on;
    subplot(1,2,2);
    plot(L, MJ_means, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Male Population by Season Length');
    xlabel('Season Length (weeks)');
    ylabel('Mean Male+Jake Density');
    grid on;
    saveas(fig, fullfile(output_dir, sprintf('pop_response_%s.png', run_id)));
    close(fig);

end

function count = actionplot(Xopt, F, ES, centerx)
    xlim([-0.075, 6.5]);
    ylim([-0.075, 6.5]);
    scatter(Xopt(:,2) + Xopt(:,3), Xopt(:,4), 50, Xopt(:,1), 'filled');
    colormap('parula');
    colorbar;
    xlabel('Total Male Density', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('Adult Female Density', 'Interpreter', 'latex', 'FontSize', 18);
    set(gca, 'FontSize', 16);
    if length(ES) >= 4
        hold on;
        plot(ES(2)+ES(3), ES(4), 'r*', 'MarkerSize', 15, 'LineWidth', 2);
        hold off;
    end
    axis square;
    count = 0;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dynamic Model Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model, results, mm, ES, aa, pp, F, J, M, svals, D, osig, omu] = ...
        PennTurkeyModel(uw_fixed, Fbar, Pbar, slope, osig, omu, use_w, cvw)

    %% addfactor tables
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
    pp    = [];
    aa    = [];

    %% Population Model Parameters
    delta  = 0.99;
    gammam = 0.409;
    gammaj = 0.653;
    gammaf = 0.542;
    cvv    = 0.0736;

    eta = 1 + 2 * Fbar * slope / Pbar;

    %% Decision Variables
    L = (0:3)';

    %% Oak Mast: build truncated distribution then clip to state grid
    O_full  = (0:26)';
    po_full = normpdf(O_full, omu, osig);
    po_full(O_full <= 5)  = sum(po_full(O_full <= 5));
    po_full(O_full >= 15) = sum(po_full(O_full >= 15));
    po = po_full(6:16) / sum(po_full(6:16));
    O  = O_full(6:16);

    %% Harvest Rate: conditioned on observed mast AND season length
    [L_grid, O_grid] = ndgrid(0:3, O);
    H_vals   = max(0.01, 0.07 - 0.00686 * O_grid + 0.0175 * L_grid);
    H_lookup = H_vals(:);

    H_cpd.type       = 'd';
    H_cpd.parameters = [];
    H_cpd.cpt        = H_vals(:)';
    H_cpd.values     = 1;
    H_cpd.q          = [];
    H_cpd.ztype      = 'u';
    H_cpd.simfunc    = @(u, idx) H_lookup(idx);

    %% Mast RV for i.i.d. transition
    o_rv.type       = 'd';
    o_rv.parameters = [];
    o_rv.values     = O(:);
    o_rv.q          = po(:);
    o_rv.cpt        = cumsum(po(:));

    %% Poult production - no mast effect
    if use_w
        Ps = @(F_val, w) 2 * eta * Pbar ./ ...
            (eta + 1 + (eta - 1) * (F_val / Fbar).^eta) .* w;
    else
        Ps = @(F_val, w) 2 * eta * Pbar ./ ...
            (eta + 1 + (eta - 1) * (F_val / Fbar).^eta);
    end

    %% Remaining Transition Functions
    Ms = @(Md, v)             gammam * Md .* v;
    Js = @(Jd, v)             gammaj * Jd .* v;
    Fs = @(Fd, v)             gammaf * Fd .* v;
    Md = @(Ms_val, Js_val)    Ms_val + Js_val;
    Jd = @(Fs_val, Pa, H_val) (1 - H_val) .* Fs_val .* Pa / 2;
    Fd = @(Fs_val, Pa, H_val) (1 - H_val) .* Fs_val .* (1 + Pa / 2);

    %% Random Variables
    vparams = Burr3mom2param(1, cvv, 0);
    v = rvdef('burr3', vparams, 25);

    if cvw > 0
        wparams = Burr3mom2param(1, cvw, 0);
        w = rvdef('burr3', wparams, 25);
    else
        wparams = Burr3mom2param(1, 1e-6, 0);
        w = rvdef('burr3', wparams, 25);
    end

    %% State Variable Grids
    minpop = 0.00;
    M = linspace(minpop, 4, 41)';
    J = linspace(minpop, 4, 41)';
    F = linspace(minpop, 7, 71)';

    %% Utility Function
    ut = @(L_val, Md_val, Jd_val, Fd_val) arrayfun(@(l, md, jd, fd) ...
        ifthenelse(fd == F(1) && l > 0, -inf, ...
        ((l + addfactor1_table(l+1)).^uw_fixed .* ...
         (md + jd + addfactor2_table(l+1)).^(1-uw_fixed)) .* (fd > F(1)) - l * 1e-10), ...
        L_val, Md_val, Jd_val, Fd_val);

    %% Decision Diagram
    D = [];
    D = add2diagram(D, 'L',       'a', true, {},              L,       [0.1278, 0.3713], []);
    D = add2diagram(D, 'Mj',      's', true, {},              M,       [0.1196, 0.7941], []);
    D = add2diagram(D, 'Jj',      's', true, {},              J,       [0.1125, 0.7143], []);
    D = add2diagram(D, 'Fj',      's', true, {},              F,       [0.1172, 0.6160], []);
    D = add2diagram(D, 'Oj',      's', true, {},              O,       [0.1200, 0.5000], []);
    D = add2diagram(D, 'v',       'c', true, {},              v,       [0.3005, 0.9365], []);
    D = add2diagram(D, 'w',       'c', true, {},              w,       [0.3282, 0.4451], []);
    D = add2diagram(D, 'H',       'c', true, {'L','Oj'},      H_cpd,   [0.5181, 0.3684], [5 1; 5 1]);
    D = add2diagram(D, 'Onext',   'c', true, {},              o_rv,    [0.6500, 0.5000], []);
    D = add2diagram(D, 'Oj+',     'f', true, {'Onext'},       @(x) x, [0.7700, 0.5000], [5 1]);
    D = add2diagram(D, 'Ms',      'c', true, {'Mj','v'},      Ms,      [0.4804, 0.7865], [5 1; 5 1]);
    D = add2diagram(D, 'Js',      'c', true, {'Jj','v'},      Js,      [0.4810, 0.7025], [5 1; 5 1]);
    D = add2diagram(D, 'Fs',      'c', true, {'Fj','v'},      Fs,      [0.4751, 0.6156], [5 1; 5 1]);
    D = add2diagram(D, 'Ps',      'c', true, {'Fj','w'},      Ps,      [0.4760, 0.4859], [5 1; 5 1]);
    D = add2diagram(D, 'Mj+',     'f', true, {'Ms','Js'},     Md,      [0.7376, 0.7955], [5 1; 5 1]);
    D = add2diagram(D, 'Jj+',     'f', true, {'Fs','Ps','H'}, Jd,      [0.7477, 0.7058], [5 1; 5 1; 5 1]);
    D = add2diagram(D, 'Fj+',     'f', true, {'Fs','Ps','H'}, Fd,      [0.7706, 0.6217], [5 1; 5 1; 5 1]);
    D = add2diagram(D, 'utility', 'u', true, {'L','Mj','Jj','Fj'}, ut, [0.3649, 0.1564], [5 1; 5 1; 5 1; 5 1]);

    %% Model and MDP Options
    doptions = struct('d', delta, 'cleanup', 2, 'reps', 0, 'chunk', 1000, 'print', 1, 'ptype', 0);
    model = d2model(D, doptions);

    moptions = struct('print', 0);
    results = mdpsolve(model, moptions);

    if isfield(results, 'Ixopt')
        results.Xopt = model.X(results.Ixopt, :);
    else
        results.Xopt = [];
    end

    %% Attempt to compute pstar from model.P if mdpsolve did not return it
    if isempty(results.pstar)
        if ~isempty(model.P) && ~isa(model.P, 'function_handle')
            try
                results.pstar = model.P(:, results.Ixopt)';
                fprintf('pstar computed from model.P: %s\n', mat2str(size(results.pstar)));
            catch
                results.pstar = [];
            end
        else
            fprintf('model.P is EV function - using simulation for mm and aa\n');
        end
    end

    %% Post-solution Processing
    if ~isempty(results.pstar)
        ns = size(results.pstar, 1);
        pp = ones(ns, 1) / ns;
        for i = 1:500
            pp = results.pstar * pp;
        end
        mm = marginals(pp, [length(M) length(J) length(F) length(O)]);
        ES = pp' * [results.Xopt Ps(results.Xopt(:,4), 1)];
        if nargout >= 5
            aa = zeros(1, 4);
            for i = 0:3
                aa(i+1) = sum(pp(results.Xopt(:,1) == i));
            end
        end
    else
        % Simulation fallback for mm, aa, ES
        fprintf('Using simulation fallback for mm and aa...\n');
        reps_sim = 100000;
        % O_init = O(round(length(O)/2));
        O_init = O(max(1, min(length(O), round((omu - O(1)) / (O(2) - O(1))) + 1)));
        SS = dsim(D, ones(reps_sim,1) * [1.5, 1.5, 3, O_init], 200, ...
                  model.X(results.Ixopt, 1), [], [], 0);
        svals = linspace(0, max(F), 501)';
        mm = cell(1, 4);
        [mm{1}, ~] = ksdensity(SS{2}(:,end), svals);
        [mm{2}, ~] = ksdensity(SS{3}(:,end), svals);
        [mm{3}, ~] = ksdensity(SS{4}(:,end), svals);
        mm{4} = po;
        aa = zeros(1, 4);
        for i = 0:3
            aa(i+1) = mean(SS{1}(:,end) == i);
        end
        ES = [mean(SS{1}(:,end)) mean(SS{2}(:,end)) mean(SS{3}(:,end)) ...
              mean(SS{2}(:,end)+SS{3}(:,end)) mean(SS{4}(:,end))];
    end

end