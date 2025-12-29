%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%X
%% A Decision Framework for Balancing Hunting Opportunity 
% and Population Abundance in Wild Turkey Management
% Winter et al. 20XX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%X
% Decision made at (January, April, September) time-point
% January = mean mast and 'guess' at recruitment
% April = predicted recruitment metrics based on weather.
% September = known recruitment and mast 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Automated Decision Framework Analysis - All Combinations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

% Add MDPSolve path
addpath(genpath('C:/Users/vaw5154/OneDrive - The Pennsylvania State University/PhD/PSUTurkey/turkey_SDP/'))

% **USER OPTIONS**
save_plots = false;  % Set to false to skip saving plots
save_results = true;  % Set to false to skip saving utility results

% Create output directory
output_dir = 'Results/Utility_Results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **USER OPTIONS**
%% Define all combinations to test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
months = {'January', 'April', 'September'};
trends = {'Increase', 'Stable', 'Decrease'};
scenarios = {'A', 'B', 'C'};
weathers = {'warm', 'cold'};  % Only used for April

% Fixed parameters
omu = 12.998; % Unobserved mast
% omu = 20; % High mast [Sept only]
% omu = 5;  % Low mast [Sept only]

slope = 0.2;
uw_fixed = 0.2;
L = (0:3)';

% Simulation parameters
reps = 100000;
time_steps = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main loop through all combinations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_results
    all_results = struct();
    result_counter = 1;
end

for m = 1:length(months)
    scenario_month = months{m};
    
    % Weather loop (only matters for April)
    if strcmp(scenario_month, 'April')
        weather_list = weathers;
    else
        weather_list = {''};
    end
    
    for w = 1:length(weather_list)
        april_weather = weather_list{w};
        
        for t = 1:length(trends)
            scenario_trend = trends{t};
            
            for s = 1:length(scenarios)
                scenario_0 = scenarios{s};
                
                % Create unique identifier
                if strcmp(scenario_month, 'April')
                    run_id = sprintf('%s_%s_%s_%s', scenario_month, april_weather, scenario_trend, scenario_0);
                else
                    run_id = sprintf('%s_%s_%s', scenario_month, scenario_trend, scenario_0);
                end
                
                fprintf('\n========================================\n');
                fprintf('Running: %s\n', run_id);
                fprintf('========================================\n');
                
                try
                    % Set parameters based on month and trend
                    [Fbar, Pbar, osig, use_w] = setParameters(scenario_month, scenario_trend, april_weather);

                    % Load utility weights
                    fname = sprintf('Data/norm_weights_popsize_scenario_%s.csv', scenario_0);
                    T = readtable(fname, 'VariableNamingRule', 'preserve');
                    row = T(strcmp(T.("Population.Size"), scenario_trend), :);
                    setGlobalWeights(row.normalized_weights);
                    
                    % Run model
                    [model, results, mm, ES, ~, pp, F, J, M, ~, D, ~, ~] = ...
                        PennTurkeyModel(uw_fixed, Fbar, Pbar, slope, osig, omu, use_w);
                    
                    % Run simulations for each season length
                    [F_means, MJ_means] = runSimulations(D, L, reps, time_steps);

                    % Save utility results if requested
                    if save_results
                        all_results(result_counter).run_id = run_id;
                        all_results(result_counter).month = scenario_month;
                        all_results(result_counter).trend = scenario_trend;
                        all_results(result_counter).scenario = scenario_0;
                        all_results(result_counter).weather = april_weather;
                        all_results(result_counter).value_function = results.v;  % Full vector
                        all_results(result_counter).stationary_dist = pp;  % If available

                        % Calculate expected utility at steady-state
                        if ~isempty(ES) && length(ES) >= 4 && ~isempty(results.v)
                            % Find state closest to expected steady-state
                            distances = sqrt((model.X(:,2) - ES(2)).^2 + ...
                                           (model.X(:,3) - ES(3)).^2 + ...
                                           (model.X(:,4) - ES(4)).^2);
                            [~, closest_idx] = min(distances);
                            all_results(result_counter).expected_utility = results.v(closest_idx);
                        elseif ~isempty(pp) && ~isempty(results.v)
                            % Use long-run expected value from stationary distribution
                            all_results(result_counter).expected_utility = pp' * results.v;
                        elseif ~isempty(results.v)
                            % Fallback: use mean value
                            all_results(result_counter).expected_utility = mean(results.v);
                        else
                            warning('Could not calculate expected utility for %s', run_id);
                            all_results(result_counter).expected_utility = NaN;
                        end
                        
                        all_results(result_counter).optimal_action = mode(model.X(results.Ixopt, 1));
                        all_results(result_counter).F_means = F_means;
                        all_results(result_counter).MJ_means = MJ_means;
                        all_results(result_counter).ES = ES;
                        result_counter = result_counter + 1;
                    end
                    
                    % Save plots if requested
                    if save_plots
                        savePlots(results, mm, ES, F, L, F_means, MJ_means, run_id, output_dir);
                    end
                    
                catch ME
                    warning('Error in run %s: %s', run_id, ME.message);
                    fprintf('Error details: %s\n', ME.getReport());
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save all results to file if requested
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_results && exist('all_results', 'var') && ~isempty(all_results)
    % Save full structure as MAT file (keeps arrays)
    save(fullfile(output_dir, 'all_utility_results.mat'), 'all_results');
    
    % Create simplified structure for CSV (scalars only)
    results_simple = struct();
    for i = 1:length(all_results)
        results_simple(i).run_id = all_results(i).run_id;
        results_simple(i).month = all_results(i).month;
        results_simple(i).trend = all_results(i).trend;
        results_simple(i).scenario = all_results(i).scenario;
        results_simple(i).weather = all_results(i).weather;
        results_simple(i).expected_utility = all_results(i).expected_utility;
        results_simple(i).optimal_action = all_results(i).optimal_action;
    end
    
    % Convert to table and save
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

function [Fbar, Pbar, osig, use_w] = setParameters(scenario_month, scenario_trend, april_weather)
    % Set Fbar and Pbar based on trend
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
    
    % Adjust Pbar for April based on weather
    if strcmp(scenario_month, 'April')
        if strcmp(april_weather, 'cold')
            Pbar = Pbar_base - 0.035;
        else  % warm
            Pbar = Pbar_base + 0.039;
        end
    else
        Pbar = Pbar_base;
    end
    
    % Set osig and use_w based on month
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

function savePlots(results, mm, ES, F, L, F_means, MJ_means, run_id, output_dir)
    % 1. Action plot
    if ~isempty(results.Xopt)
        fig = figure('Visible', 'off');
        try
            actionplot(results.Xopt, F, ES, 1);
            %title(strrep(run_id, '_', '-'), 'interpreter', 'latex');
            grid on;
            saveas(fig, fullfile(output_dir, sprintf('action_%s.png', run_id)));
        catch
            warning('Could not create action plot for %s', run_id);
        end
        close(fig);
    end
    
    % 2. Season length distribution
    if ~isempty(mm) && ~isempty(mm{1})
        fig = figure('Visible', 'off');
        mm_trimmed = mm{1}(1:4);
        lrharvplot_fixed(mm_trimmed);
       % title(strrep(run_id, '_', '-'), 'interpreter', 'latex');
        saveas(fig, fullfile(output_dir, sprintf('season_dist_%s.png', run_id)));
        close(fig);
    end
    
    % 3. Population response
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

function lrharvplot_fixed(EL, width, color)
    if nargin < 2, width = 2; end
    if nargin < 3, color = [0 0 0]; end
    
    L_plot = (0:length(EL)-1)';
    plot(L_plot, EL, '.-', 'LineWidth', width, 'Color', color, 'MarkerSize', 20);
    ylim([0 max([1, max(EL)*1.1])]);
    xlim([-0.5, length(L_plot)-0.5]);
    set(gca, 'XTick', 0:1:length(L_plot), 'FontSize', 16);

    xlabel('Season Length (weeks)', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('Probability', 'Interpreter', 'latex', 'FontSize', 18);

    grid on;
end

function count = actionplot(Xopt, F, ES, centerx)
    xlim([-0.075,6.5]);
    ylim([-0.075,6.5]);
    
    scatter(Xopt(:,2) + Xopt(:,3), Xopt(:,4), 50, Xopt(:,1), 'filled');
    colormap('parula');
    colorbar;
    cb.FontSize = 16;

    xlabel('Total Male Density', 'interpreter', 'latex', 'FontSize', 18);
    ylabel('Adult Female Density', 'interpreter', 'latex', 'FontSize', 18);
    
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

% Dynamic model function
function [model, results, mm, ES, aa, pp, F, J, M, svals, D, osig, omu] = PennTurkeyModel(uw_fixed, Fbar, Pbar, slope, osig, omu, use_w)
      %% Parameter Definitions and Adjustments
    addfactor = 0.1;
    
    % Calculate action-specific addfactor1 values using hunter preferences from getUtilityWeight
    addfactor1_table = zeros(4,1);
    addfactor2_table = zeros(4,1);
    for i = 0:3
        hunter_pref = getUtilityWeight(i);  % Get hunter preference weight for this action
        
        % Calculate addfactor1 based on hunter preferences (ensures U(1,y)=2U(0,y))
        if hunter_pref == 1
            addfactor1_table(i+1) = 0.1;  % Avoid division by zero
        else
            addfactor1_table(i+1) = 1/(2^(1/hunter_pref)-1);
        end
        
        % Calculate addfactor2 based on hunter preferences (ensures U(x,1)=2U(x,0))
        if hunter_pref == 0
            addfactor2_table(i+1) = 0.1;  % Avoid division by zero
        else
            addfactor2_table(i+1) = 1/(2^(1/(1-hunter_pref))-1);
        end
    end
    
    % Reset addfactor2 to zero if desired (as in Fackler original code)
    addfactor2_table = zeros(4,1);  % Reset all to zero
    
    % Initialize state grids and auxiliary variables
    svals = [];
    pp    = [];
    aa    = [];
    
    %% Discount Factor and Population Model Parameters
    delta   = 0.99;      % Discount factor
    
    % Summer survival parameters
    gammam  = 0.409;     % Male spring/summer mean survival 
    gammaj  = 0.653;     % Jake spring/summer mean survival 
    gammaf  = 0.542;     % Female spring/summer mean survival
    
    % Coefficients of variation for survival rates
    cvv     = 0.0736;         % Summer survival variation
    cvw     = 1.116 / 4.96;   % Variation for poults/hens survival

    % Fecundity shape parameter (derived from Fbar, slope, and Pbar)
    eta = 1 + 2 * Fbar * slope / Pbar;

    %% Decision Variables and Harvest Rate Computation
    % Season length decision variable (in weeks)
    L = (0:3)';
    
    % Oak mast levels
    O = (0:26)';
    
    % Probability distribution of oak mast levels (assumed normal)
    po = normpdf(O, omu, osig);
    
    % Truncate distribution: Sum probabilities for O<=5 and O>=15 
    po(O <= 5) = sum(po(O <= 5));
    po(O >= 15) = sum(po(O >= 15));
    po = po(6:16) / sum(po(6:16));
    
    % Harvest rate table: a minimum of 0.01 is enforced 
    H = [zeros(11,1) max(0.01,(0.07 - 0.00686*O(6:16) + 0.0175*(1:4)))];
    eH = H' * po;  % Expected harvest over mast levels
    
    %% Transition Functions Definitions
    % Summer state transitions:
    Ms = @(Md, v) gammam * Md .* v;     % September male population
    Js = @(Jd, v) gammaj * Jd .* v;     % September jake population
    Fs = @(Fd, v) gammaf * Fd .* v;     % September female population
    H  = @(L)    eH(L+1);               % harvest rate (from table)

    % Poult production function with density-dependent fecundity
    % Ps = @(F, w) 2 * eta * Pbar ./ (eta + 1 + (eta - 1) * (F / Fbar).^eta) .* w; % 1 for september poults per hen
    % Poult production function
    if use_w
        Ps = @(F, w) 2 * eta * Pbar ./ (eta + 1 + (eta - 1) * (F / Fbar).^eta) .* w;
    else
        Ps = @(F, w) 2 * eta * Pbar ./ (eta + 1 + (eta - 1) * (F / Fbar).^eta) * 1;
    end
    % Winter transitions:
    Md = @(Ms, Js) Ms + Js;                               % December male population
    Jd = @(Fs, Pa, H) (1 - H) .* Fs .* Pa / 2;            % December jake population
    Fd = @(Fs, Pa, H) (1 - H) .* Fs .* (1 + Pa / 2);      % December female population
     
    %% Random Noise Variables for Survival and Reproduction
    vparams = Burr3mom2param(1, cvv, 0);
    wparams = Burr3mom2param(1, cvw, 0);
    v = rvdef('burr3', vparams, 25);            % Summer survival
     w = rvdef('burr3', wparams, 25);            % Poults per hen
    
    %% Discretization of the State Variables
    minpop = 0.00;
    % M = linspace(minpop, 1.5, 41)';
    % J = linspace(minpop, 1.0, 41)';
    % F = linspace(minpop, 3.0, 71)';

     M = linspace(minpop, 4, 41)';  % Male population grid
     J = linspace(minpop, 4, 41)';  % Jake population grid
     F = linspace(minpop, 7, 71)';  % Female population grid
      
       %% Integrated Utility Function with Hunter Preference-Influenced Scaling
    % Use fixed uw_fixed for consistent trade-off, but hunter preferences influence scaling
    utility = @(L, Md, Jd, Fd) arrayfun(@(l, md, jd, fd) ...
        ((l + addfactor1_table(l+1)).^uw_fixed .* ...
         (md + jd + addfactor2_table(l+1)).^(1-uw_fixed)) .* (fd > F(1)) - l * 1e-10, ...
        L, Md, Jd, Fd);

    % Handle minimum female population case  
    ut = @(L, Md, Jd, Fd) arrayfun(@(l, md, jd, fd) ...
        ifthenelse(fd == F(1) && l > 0, -inf, ...
        ((l + addfactor1_table(l+1)).^uw_fixed .* ...
         (md + jd + addfactor2_table(l+1)).^(1-uw_fixed)) .* (fd > F(1)) - l * 1e-10), ...
        L, Md, Jd, Fd);
     
     
    %% Building the Decision Diagram
    D = [];
    % January
    D = add2diagram(D, 'L', 'a', true, {}, L, [0.1278, 0.3713], []);
    D = add2diagram(D, 'Mj', 's', true, {}, M, [0.1196, 0.7941], []);
    D = add2diagram(D, 'Jj', 's', true, {}, J, [0.1125, 0.7143], []);
    D = add2diagram(D, 'Fj', 's', true, {}, F, [0.1172, 0.6160], []);
    % September t 
    D = add2diagram(D, 'v', 'c', true, {}, v, [0.3005, 0.9365], []);
    D = add2diagram(D, 'w', 'c', true, {}, w, [0.3282, 0.4451], []);
    D = add2diagram(D, 'H', 'c', true, {'L'}, H, [0.5181, 0.3684], [5 1]);
    D = add2diagram(D, 'Ms', 'c', true, {'Mj','v'}, Ms, [0.4804, 0.7865], [5 1; 5 1]);
    D = add2diagram(D, 'Js', 'c', true, {'Jj','v'}, Js, [0.4810, 0.7025], [5 1; 5 1]);
    D = add2diagram(D, 'Fs', 'c', true, {'Fj','v'}, Fs, [0.4751, 0.6156], [5 1; 5 1]);
    D = add2diagram(D, 'Ps', 'c', true, {'Fj','w'}, Ps, [0.4760, 0.4859], [5 1; 5 1]);
    % December t + 1 
    D = add2diagram(D, 'Mj+', 'f', true, {'Ms','Js'}, Md, [0.7376, 0.7955], [5 1; 5 1]);
    D = add2diagram(D, 'Jj+', 'f', true, {'Fs','Ps','H'}, Jd, [0.7477, 0.7058], [5 1; 5 1; 5 1]);
    D = add2diagram(D, 'Fj+', 'f', true, {'Fs','Ps','H'}, Fd, [0.7706, 0.6217], [5 1; 5 1; 5 1]);
    D = add2diagram(D, 'utility', 'u', true, {'L','Mj','Jj','Fj'}, ut, [0.3649, 0.1564], [5 1; 5 1; 5 1; 5 1]);

    %% Model and MDP Options
    doptions = struct('d', delta, 'cleanup', 2, 'reps', 0, 'chunk', 1000, 'print', 1, 'ptype', 0);
    model = d2model(D, doptions);
    
    moptions = struct('print', 0);
    results = mdpsolve(model, moptions);
    
    if isfield(results, 'Ixopt')
        Xopt = model.X(results.Ixopt, :);
        results.Xopt = Xopt;  % Store for consistency
    else
        results.Xopt = [];
    end
     
    %% Post-solution Processing: Marginals and Expected State Values
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
        % If pstar is empty, run simulation to estimate marginals and expectations
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