clc;
clear
close

%% PRCC_R0_script.m
% This script performs a global sensitivity analysis on the computed R0 by:
%   1. Sampling key model parameters using Latin Hypercube Sampling (LHS)
%   2. Computing R0 for each parameter set
%   3. Calculating Partial Rank Correlation Coefficients (PRCC)
%   4. Plotting PRCC values and contour plots of R0

%% Load Data and Define Baseline Parameters
Polymod     = load("MixingData.mat");
varyingpara = load("parameters.mat");
prob        = load("Probab.mat");
Tdata       = readtable('test_data.xlsx','Sheet','Jan_SepTest');
Pop         = load("EnglandDecilePop.mat");

N     = Pop.TotalPop;
D     = sum(N,2)';
N_age = sum(N,1);
total = sum(D);

Mix_H = Polymod.UK_from_toH;
Mix_S = Polymod.UK_from_toS;
Mix_W = Polymod.UK_from_toW;
Mix_O = Polymod.UK_from_toO;

a = varyingpara.Susceptibility;
h = varyingpara.Detection;

para = struct('tau',1, 'epsilon',2/100, 'p',0.8526, 'rho',0.5, 'theta',0.5, ...
    'delta',prob.Sympt_2_hosp, 'gamma',1/14, 'gamma_2',1/14, 'd',prob.Hosp_2_Death, ...
    'sigma',1/3, 'a',a, 'h',h, 'N',Polymod.UK_PP, 'beta',2.5, 'iota',0.4, 'N_total',total);

pi_val = 4.25e-6;

%% Build Age and Social Mixing Matrices
AgeMat0 = (Mix_H + Mix_S + Mix_W + Mix_O) .* (para.a .* para.h');
x = 10:-1:1;
epi = 0.3;
soc = SocioContMix(x, epi, total, D); 
soc = soc / sum(soc(:));

% %% Global Sensitivity via Latin Hypercube Sampling
 % paramNames = {'beta', 'rho', 'gamma', 'delta', 'iota', 'pi_val'};
 paramNames = {'beta', 'rho', 'gamma', 'delta', 'iota', 'pi_val', 'lambda_scale'};
 
nParams = length(paramNames);
nSim = 1000;

baseline = [para.beta, para.rho, para.gamma, 1.0, para.iota, pi_val, 1.0];
lower = baseline .* [0.8, 0.8, 0.8, 0.3, 0.8, 0.8, 0.8];
upper = baseline .* [1.2, 1.2, 1.2, 0.7, 1.2, 1.2, 1.2];

LHS = lhsdesign(nSim, nParams);
paramSamples = repmat(lower, nSim, 1) + LHS .* repmat((upper - lower), nSim, 1);

R0_samples = zeros(nSim, 1);

for i = 1:nSim
    para_sample = para;
    pi_sample = pi_val;
    lambda_scale = 1.0;

    for j = 1:nParams
        param = paramNames{j};
        val = paramSamples(i, j);
        if strcmp(param, 'pi_val')
            pi_sample = val;
        elseif strcmp(param, 'delta')
            para_sample.delta = para.delta*val;
        elseif strcmp(param, 'lambda_scale')
            lambda_scale = val;
        else
            para_sample.(param) = val;
        end
    end

    R0_samples(i) = computeR0_from_para(para_sample, pi_sample, soc, AgeMat0 * lambda_scale);
end

%% Compute PRCC
rankedParams = tiedrank(paramSamples);
rankedR0 = tiedrank(R0_samples);
prcc_values = zeros(nParams, 1);

for j = 1:nParams
    otherIdx = setdiff(1:nParams, j);
    X_other = [ones(nSim,1), rankedParams(:, otherIdx)];
    res_j = rankedParams(:, j) - X_other * regress(rankedParams(:, j), X_other);
    res_R0 = rankedR0 - X_other * regress(rankedR0, X_other);
    % [~,~,res_j] = regress(rankedParams(:, j), X_other);
    % [~,~,res_R0] = regress(rankedR0, X_other);
    prcc_values(j) = corr(res_j, res_R0);
end

%% Plot PRCC Results
% Define Greek/LaTeX labels for parameters, including 'pi_val' and 'iota'
greek_labels_map = containers.Map(...
    {'beta', 'epsilon', 'rho', 'theta', 'gamma', 'sigma', 'iota', ...
     'delta', 'd', 'a', 'h', 'pi_val', 'lambda_scale'}, ...
    {'$\hat{q}$', '$\epsilon$', '$\rho$', '$\theta$', '$\gamma$', ...
     '$\sigma$', '$\iota$', '$\delta$', '$d$', '$a$', '$h$', '$\pi_{i}$', '$\lambda$'});

xtick_labels = cell(size(paramNames));
for i = 1:length(paramNames)
    name = paramNames{i};
    if isKey(greek_labels_map, name)
        xtick_labels{i} = greek_labels_map(name);
    else
        xtick_labels{i} = ['$', strrep(name, '_', '\_'), '$'];
    end
end

figure;
bar(prcc_values, 'FaceColor', [0.3 0.7 0.9]);
ylim([-1, 1]);
grid on;
xlabel("Parameters", 'Interpreter', 'latex');
ylabel("Partial Rank Correlation", 'Interpreter', 'latex');

set(gca, 'XTick', 1:length(paramNames), ...
         'XTickLabel', xtick_labels, ...
         'TickLabelInterpreter', 'latex', ...
         'FontSize', 18);


% Contour Plots for Selected Pairs
paramNames = {'beta', 'rho', 'gamma', 'delta'};

greek_labels_map = containers.Map(...
    {'beta', 'epsilon', 'rho', 'theta', 'gamma', 'sigma', 'iota', ...
     'delta', 'd', 'a', 'h', 'pi_val', 'lambda_scale'}, ...
    {'$\hat{q}$', '$\epsilon$', '$\rho$', '$\theta$', '$\gamma$', ...
     '$\sigma$', '$\iota$', '$\delta$', '$d$', '$a$', '$h$', '$\pi_{i}$', '$\lambda$'});

xtick_labels = cell(size(paramNames));
for i = 1:length(paramNames)
    name = paramNames{i};
    if isKey(greek_labels_map, name)
        xtick_labels{i} = greek_labels_map(name);
    else
        xtick_labels{i} = ['$', strrep(name, '_', '\_'), '$'];
    end
end

nGrid = 50; scaling = 0.5;
pairs = nchoosek(1:length(paramNames), 2);
nPairs = size(pairs, 1);
nCols = 3; nRows = ceil(nPairs / nCols);
figure('Units','normalized','Position',[0.05 0.05 0.9 0.85]);
tiledlayout(nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');
labels = arrayfun(@(x) char('A'+x-1), 1:nPairs, 'UniformOutput', false);

for k = 1:nPairs
    idx1 = pairs(k,1); idx2 = pairs(k,2);
    p1_name = paramNames{idx1}; p2_name = paramNames{idx2};
    p1_base = para.(p1_name); p2_base = para.(p2_name);

    p1_vals = strcmp(p1_name,'delta') * linspace(0.2,1.0,nGrid) + ~strcmp(p1_name,'delta') * linspace(p1_base*(1-scaling), p1_base*(1+scaling), nGrid);
    % p2_vals = strcmp(p2_name,'delta') * linspace(0.2,1.0,nGrid) + ~strcmp(p2_name,'delta') * linspace(p2_base*(1-scaling), p2_base*(1+scaling), nGrid);

    if strcmp(p2_name, 'delta')
        p2_vals = linspace(0.2, 1.0, nGrid);
    else
        p2_vals = linspace(p2_base * (1 - scaling), p2_base * (1 + scaling), nGrid);
    end

    % if strcmp(p1_name, 'gamma')
    %     p1_vals = linspace(1/21, 1/7, nGrid);
    % % s * (1 - scaling), p1_base * (1 + scaling), nGrid);
    % end
    % 
    %  if strcmp(p2_name, 'gamma')
    %     p2_vals = linspace(1/21, 1/7, nGrid);
    %  end


    % if strcmp(p2_name, 'gamma')
    %     p2_vals = linspace(1/14, 1/7, nGrid);
    % else
    %     p2_vals = linspace(p2_base * (1 - scaling), p2_base * (1 + scaling), nGrid);
    % end



    R0_grid = zeros(nGrid, nGrid);
    for i = 1:nGrid
        for j = 1:nGrid
            temp_para = para;
            temp_para.(p1_name) = strcmp(p1_name,'delta') * (para.delta * p1_vals(i)) + ~strcmp(p1_name,'delta') * p1_vals(i);
            temp_para.(p2_name) = strcmp(p2_name,'delta') * (para.delta * p2_vals(j)) + ~strcmp(p2_name,'delta') * p2_vals(j);
            R0_grid(j,i) = computeR0_from_para(temp_para, pi_val, soc, AgeMat0);
        end
    end

    ax = nexttile;
    [X, Y] = meshgrid(p1_vals, p2_vals);
    contourf(X, Y, R0_grid, 20, 'LineColor', 'none');
    colormap(parula); colorbar;


     if strcmp(p1_name, 'beta')
         xlabel('$\hat{q}$', 'Interpreter', 'latex');
        elseif isKey(greek_labels_map, p1_name)
            xlabel(['$', greek_labels_map(p1_name), '$'], 'Interpreter', 'latex');
        else
            xlabel(strrep(p1_name, '_', '\_'), 'Interpreter', 'latex');
     end

        if strcmp(p2_name, 'beta')
            ylabel('$\hat{q}$', 'Interpreter', 'latex');
        elseif isKey(greek_labels_map, p2_name)
            ylabel(['$', greek_labels_map(p2_name), '$'], 'Interpreter', 'latex');
        else
            ylabel(strrep(p2_name, '_', '\_'), 'Interpreter', 'latex');
        end


    % p2_label = strcmp(p2_name, 'beta') * '\hat{q}' + ~strcmp(p2_name, 'beta') * strrep(p2_name, '_', '\');

    % xlabel(['$', p1_label, '$'], 'Interpreter','latex','FontSize',18);
    % ylabel(['$', p2_label, '$'], 'Interpreter','latex','FontSize',18);
    text(0.02, 0.95, labels{k}, 'Units','normalized', 'FontWeight','bold', 'FontSize',18, ...
        'BackgroundColor','w', 'Margin',1, 'EdgeColor','k', 'LineWidth',1.5, 'Parent',ax);
end

sgtitle('Contour Plots of $R_0$ vs Selected Parameter Pairs', ...
    'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'latex');
print(gcf, 'R0_Contour_SelectedPairs', '-dpng', '-r300');

% 

%% Function to Compute R0
function R0 = computeR0_from_para(para, pi_val, soc, AgeMat0)

    
    term1 = ((1 - pi_val) * para.rho) / para.gamma;
    term2 = ((1 - pi_val) * (1 - para.rho)) ./ (para.delta + para.gamma);
    term3 = (para.iota * pi_val) ./ (para.delta + para.gamma);
    scalar_sum = term1 + term2 + term3;

    n_decile = size(soc, 1);
    n_age = size(AgeMat0, 1);
    AS_0 = kron(soc, AgeMat0);

    riskgroups = zeros(n_decile,n_decile);
    for i = 1:n_decile
        for j = 1:n_decile
            rows = ((i-1)*n_age + 1):(i*n_age);
            cols = ((j-1)*n_age + 1):(j*n_age);
            % disp(size(AS_0(rows, cols)))
            block = para.beta .* (AS_0(rows, cols) .* (scalar_sum));
            riskgroups(i,j) = max(real(eig(block)));
        end
    end
    
    R0 = max(eig(riskgroups));
end


% 
% % %% Contour plot two
% 
% % Contour Plots for Selected Pairs (Including lambda_scale)
% paramNames = {'lambda_scale','rho', 'gamma', 'delta'};
% nGrid = 50; scaling = 0.2;
% pairs = nchoosek(1:length(paramNames), 2);
% nPairs = size(pairs, 1);
% nCols = 3; nRows = ceil(nPairs / nCols);
% 
% % Define baseline for lambda_scale
% lambda_base = 1.0;
% 
% figure('Units','normalized','Position',[0.05 0.05 0.9 0.85]);
% tiledlayout(nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');
% labels = arrayfun(@(x) char('A'+x-1), 1:nPairs, 'UniformOutput', false);
% 
% for k = 1:nPairs
%     idx1 = pairs(k,1); idx2 = pairs(k,2);
%     p1_name = paramNames{idx1}; p2_name = paramNames{idx2};
% 
%     % Set parameter ranges
%     if strcmp(p1_name, 'lambda_scale')
%         p1_vals = linspace(0.1, 1.0, nGrid);
%     elseif strcmp(p1_name, 'gamma')
%         p1_vals = linspace(1/21, 1/7, nGrid);
%     elseif strcmp(p1_name, 'delta')
%         p1_vals = linspace(0.2, 1.0, nGrid);
%     else
%         p1_vals = linspace(para.(p1_name)*(1-scaling), para.(p1_name)*(1+scaling), nGrid);
%     end
% 
%     if strcmp(p2_name, 'lambda_scale')
%         p2_vals = linspace(0.1, 1.0, nGrid);
%     elseif strcmp(p2_name, 'gamma')
%         p2_vals = linspace(1/21, 1/7, nGrid);
%     elseif strcmp(p2_name, 'delta')
%         p2_vals = linspace(0.2, 1.0, nGrid);
%     else
%         p2_vals = linspace(para.(p2_name)*(1-scaling), para.(p2_name)*(1+scaling), nGrid);
%     end
% 
%     R0_grid = zeros(nGrid, nGrid);
%     for i = 1:nGrid
%         for j = 1:nGrid
%             temp_para = para;
%             lambda_eff = lambda_base; % Default to baseline
%             if strcmp(p1_name, 'lambda_scale')
%                 lambda_eff = p1_vals(i);
%             elseif strcmp(p2_name, 'lambda_scale')
%                 lambda_eff = p2_vals(j);
%             else
%                 % Update parameters only if not lambda_scale
%                 if strcmp(p1_name, 'delta')
%                     temp_para.delta = para.delta * p1_vals(i);
%                 else
%                     temp_para.(p1_name) = p1_vals(i);
%                 end
%                 if strcmp(p2_name, 'delta')
%                     temp_para.delta = para.delta * p2_vals(j);
%                 else
%                     temp_para.(p2_name) = p2_vals(j);
%                 end
%             end
%             AgeMat_temp = AgeMat0 * lambda_eff; % Temporary matrix
%             R0_grid(j,i) = computeR0_from_para(temp_para, pi_val, soc, AgeMat_temp);
%         end
%     end
% 
%     ax = nexttile;
%     [X, Y] = meshgrid(p1_vals, p2_vals);
%     contourf(X, Y, R0_grid, 20, 'LineColor', 'none');
%     colormap(parula); colorbar;
% 
%     if strcmp(p1_name, 'lambda_scale')
%         xlabel('$\lambda$', 'Interpreter', 'latex');
%     elseif isKey(greek_labels_map, p1_name)
%         xlabel(['$', greek_labels_map(p1_name), '$'], 'Interpreter', 'latex');
%     else
%         xlabel(strrep(p1_name, '_', '\_'), 'Interpreter', 'latex');
%     end
% 
%     if strcmp(p2_name, 'lambda_scale')
%         ylabel('$\lambda$', 'Interpreter', 'latex');
%     elseif isKey(greek_labels_map, p2_name)
%         ylabel(['$', greek_labels_map(p2_name), '$'], 'Interpreter', 'latex');
%     else
%         ylabel(strrep(p2_name, '_', '\_'), 'Interpreter', 'latex');
%     end
% 
%     text(0.02, 0.95, labels{k}, 'Units','normalized', 'FontWeight','bold', 'FontSize',18, ...
%         'BackgroundColor','w', 'Margin',1, 'EdgeColor','k', 'LineWidth',1.5, 'Parent',ax);
% end
% 
% sgtitle('Contour Plots of $R_0$ vs Selected Parameter Pairs (using mixing matrix)', ...
%     'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'latex');
% print(gcf, 'R0_Contour_SelectedPairs_Lambda', '-dpng', '-r300');
