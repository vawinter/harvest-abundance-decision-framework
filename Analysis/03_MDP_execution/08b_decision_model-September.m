%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A Decision Framework for Balancing Hunting Opportunity and Population Abundance:
% A Case Study for Wild Turkey Management
% Winter et al. 20XX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decision made at September time-point
% September = known recruitment and mast
% THIS ONE WORKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

% Add MDPSolve path
addpath(genpath('C:/Users/vaw5154/OneDrive - The Pennsylvania State University/PhD/PSUTurkey/turkey_SDP/'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER OPTIONS
save_plots   = true;
save_results = true;
output_dir   = 'Results/Sept_MastState_Analysis';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define all combinations to test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
months    = {'September'};
trends    = {'Increase', 'Stable', 'Decrease'};
scenarios = {'A', 'B', 'C'};
weathers  = {'warm', 'cold'};
model_var = {'Low', 'High'};

% Fixed parameters
% Average mast
omu = 12.998;

% Low mast  
% omu = 5;

% High mast
% omu = 20;

slope    = 0.2;
uw_fixed = 0.2;
L        = (0:3)';

% Simulation parameters (used only if pstar is empty)
reps       = 100000;
time_steps = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_results
    all_results    = struct();
    result_counter = 1;
end


for m = 1:length(months)
    scenario_month = months{m};

    % Weather / variation lists depend on month
    if strcmp(scenario_month, 'April')
        weather_list = weathers;
        var_list     = model_var;
    elseif strcmp(scenario_month, 'September')
        weather_list = weathers;
        var_list     = {''};
    else  % January
        weather_list = {''};
        var_list     = {''};
    end

    for v = 1:length(var_list)
        april_var = var_list{v};

        for w = 1:length(weather_list)
            april_weather = weather_list{w};

            for t = 1:length(trends)
                scenario_trend = trends{t};

                for s = 1:length(scenarios)
                    scenario_0 = scenarios{s};

                    % Unique run identifier
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
                        %% Set parameters
                        [Fbar, Pbar, osig, cvw] = setParameters(scenario_month, scenario_trend, april_weather, april_var);

                        %% Load utility weights
                        fname = sprintf('Data/norm_weights_popsize_scenario_%s.csv', scenario_0);
                        T   = readtable(fname, 'VariableNamingRule', 'preserve');
                        row = T(strcmp(T.("Population.Size"), scenario_trend), :);
                        setGlobalWeights(row.normalized_weights);

                        use_w = (cvw > 0);

                        %% Solve model
                        [model, results, mm, ES, aa, pp, F, J, M, svals, D, ~, ~] = ...
                            PennTurkeyModel(uw_fixed, Fbar, Pbar, slope, osig, omu, use_w, cvw);

                        %% Compute population response from stationary distribution
                        [F_means, MJ_means] = computePopResponse(L, results, pp);

                        %% Store results
                        if save_results
                            all_results(result_counter).run_id        = run_id;
                            all_results(result_counter).month         = scenario_month;
                            all_results(result_counter).trend         = scenario_trend;
                            all_results(result_counter).scenario      = scenario_0;
                            all_results(result_counter).weather       = april_weather;
                            all_results(result_counter).model_var     = april_var;
                            all_results(result_counter).value_function = results.v;
                            all_results(result_counter).stationary_dist = pp;
                            all_results(result_counter).ES            = ES;
                            all_results(result_counter).F_means       = F_means;
                            all_results(result_counter).MJ_means      = MJ_means;
                            all_results(result_counter).optimal_action = mode(model.X(results.Ixopt, 1));

                            % Expected utility at steady-state
                            all_results(result_counter).expected_utility = ...
                                computeExpectedUtility(model, results, pp, ES);

                            result_counter = result_counter + 1;
                        end

                        %% Save plots
                        if save_plots
                            savePlots(results, aa, mm, ES, F, M, J, L, ...
                                      F_means, MJ_means, svals, run_id, output_dir);
                        end

                    catch ME
                        warning('Error in run %s: %s', run_id, ME.message);
                        fprintf('Error details: %s\n', ME.getReport());
                    end

                end  % scenarios
            end  % trends
        end  % weather
    end  % variation
end  % months

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_results && exist('all_results', 'var') && ~isempty(all_results)
    save(fullfile(output_dir, 'all_utility_results.mat'), 'all_results');

    % Scalar-only table for CSV
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
    writetable(struct2table(results_simple), fullfile(output_dir, 'utility_results.csv'));

    fprintf('\n========================================\n');
    fprintf('Results saved: %d runs to %s\n', length(all_results), output_dir);
    fprintf('========================================\n');
else
    fprintf('\nNo results to save.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Fbar, Pbar, osig, cvw] = setParameters(scenario_month, scenario_trend, april_weather, april_var)
    % Fbar and Pbar base values by trend
    switch scenario_trend
        case 'Increase', Fbar = 2.75; Pbar_base = 1.98;
        case 'Stable',   Fbar = 2.50; Pbar_base = 1.80;
        otherwise,       Fbar = 2.25; Pbar_base = 1.62;  % Decrease
    end

    % Adjust Pbar for weather (April and September only)
    if strcmp(scenario_month, 'April') || strcmp(scenario_month, 'September')
        switch april_weather
            case 'cold', Pbar = Pbar_base - 0.035;
            case 'warm', Pbar = Pbar_base + 0.039;
            otherwise,   Pbar = Pbar_base;
        end
    else
        Pbar = Pbar_base;
    end

    % osig and cvw by month
    switch scenario_month
        case 'September'
            osig = 1.0;
            cvw  = 0;
        case 'April'
            osig = 2.3;
            cvw  = (1.116/4.96) * (strcmp(april_var,'Low')*0.5 + strcmp(april_var,'High')*2.0);
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
function eu = computeExpectedUtility(model, results, pp, ES)
    if ~isempty(pp) && ~isempty(results.v)
        eu = pp' * results.v;
    elseif ~isempty(ES) && length(ES) >= 4 && ~isempty(results.v)
        distances = sqrt((model.X(:,2) - ES(2)).^2 + ...
                         (model.X(:,3) - ES(3)).^2 + ...
                         (model.X(:,4) - ES(4)).^2);
        [~, idx] = min(distances);
        eu = results.v(idx);
    elseif ~isempty(results.v)
        eu = mean(results.v);
    else
        eu = NaN;
        warning('Could not calculate expected utility.');
    end
end

%--------------------------------------------------------------------------
function [F_means, MJ_means] = computePopResponse(L, results, pp)
    % Compute stationary-distribution-weighted mean populations
    % conditional on each season length, from the MDP solution directly.
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

%--------------------------------------------------------------------------
function savePlots(results, aa, mm, ES, F, M, J, L, F_means, MJ_means, svals, run_id, output_dir)

    %% 1. Action plot (optimal season length by state)
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
        ylim([0, max(1, max(aa)*1.1)]);
        xlim([-0.5, length(aa)-0.5]);
        set(gca, 'XTick', 0:length(aa)-1, 'FontSize', 16);
        xlabel('Season Length (weeks)', 'Interpreter', 'latex', 'FontSize', 18);
        ylabel('Probability',           'Interpreter', 'latex', 'FontSize', 18);
        grid on;
        saveas(fig, fullfile(output_dir, sprintf('season_dist_%s.png', run_id)));
        close(fig);
    end

    %% 3. Long-run population density distributions
    if ~isempty(mm) && length(mm) >= 3
        % Determine x-axes
        if ~isempty(svals) && length(mm{1}) == length(svals)
            xM = svals; xJ = svals; xF = svals;
        else
            xM = M; xJ = J; xF = F;
        end

        fig = figure('Visible', 'off');
        set(fig, 'Position', [100, 100, 700, 500]);
        hold on;
        plot(xM, mm{1}, 'b-', 'LineWidth', 2, 'DisplayName', '$M$');
        plot(xJ, mm{2}, 'r-', 'LineWidth', 2, 'DisplayName', '$J$');
        plot(xF, mm{3}, 'y-', 'LineWidth', 2, 'DisplayName', '$F$');
        hold off;
        legend('Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 14);
        xlim([0, max(F)]);
        ylim([0, max([max(mm{1}), max(mm{2}), max(mm{3})]) * 1.15]);
        xlabel('Population Density',  'Interpreter', 'latex', 'FontSize', 16);
        ylabel('Probability Density', 'Interpreter', 'latex', 'FontSize', 16);
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

%--------------------------------------------------------------------------
function count = actionplot(Xopt, F, ES, centerx)
    scatter(Xopt(:,2)+Xopt(:,3), Xopt(:,4), 50, Xopt(:,1), 'filled');
    colormap('parula');
    cb = colorbar; cb.FontSize = 16;
    xlabel('Total Male Density',   'Interpreter', 'latex', 'FontSize', 18);
    ylabel('Adult Female Density', 'Interpreter', 'latex', 'FontSize', 18);
    xlim([-0.075, 6.5]);
    ylim([-0.075, 6.5]);
    set(gca, 'FontSize', 16);
    if length(ES) >= 4
        hold on;
        plot(ES(2)+ES(3), ES(4), 'r*', 'MarkerSize', 15, 'LineWidth', 2);
        hold off;
    end
    axis square;
    count = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dynamic Model Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model, results, mm, ES, aa, pp, F, J, M, svals, D, osig, omu] = ...
        PennTurkeyModel(uw_fixed, Fbar, Pbar, slope, osig, omu, use_w, cvw)

    %% Utility weight addfactor tables (one entry per season length 0-3)
    addfactor1_table = zeros(4,1);
    addfactor2_table = zeros(4,1);
    for i = 0:3
        hp = getUtilityWeight(i);
        addfactor1_table(i+1) = ifthenelse(hp == 1, 0.1, 1/(2^(1/hp)-1));
        addfactor2_table(i+1) = 0;  % fixed to zero as in original
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
    Ogrid = Ovec(6:16);   % mast values 5:15

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
    H  = @(Lval, Oval)     max(0.01, 0.07 - 0.00686.*Oval + 0.0175.*Lval);  % H(L,O)

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

    %% Solve MDP
    % ptype=1 ensures stationary distribution is computed
    doptions = struct('d', delta, 'cleanup', 2, 'reps', 0, 'chunk', 1000, 'print', 1, 'ptype', 1);
    model    = d2model(D, doptions);

    moptions = struct('print', 0);
    results  = mdpsolve(model, moptions);

    fprintf('pstar empty: %d, Ixopt size: %d\n', isempty(results.pstar), length(results.Ixopt));

    if isfield(results, 'Ixopt')
        results.Xopt = model.X(results.Ixopt, :);
    else
        results.Xopt = [];
    end

    %% Post-solution: stationary distribution path
    if ~isempty(results.pstar)
        ns = size(results.pstar, 1);
        pp = ones(ns,1) / ns;
        for i = 1:500
            pp = results.pstar * pp;
        end

       % Marginals over M, J, F
mm = marginals(pp, [length(M) length(J) length(F)]);

% Expected values
ES = pp' * [results.Xopt, Ps(results.Xopt(:,4), 1)];

% Season length distribution -- goes into aa only, do NOT overwrite mm{1}
aa = zeros(4, 1);
for i = 0:3
    aa(i+1) = sum(pp(results.Xopt(:,1) == i));
end
% mm{1} is left as the M marginal from marginals()

    else
        %% Fallback: simulation (should not be reached with ptype=1)
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