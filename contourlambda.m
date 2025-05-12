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

para = struct('tau',1, 'epsilon',2/100, 'p',0.8526, 'rho',0.6, 'theta',0.4, ...
    'delta',prob.Sympt_2_hosp, 'gamma',1/14, 'gamma_2',1/14, 'd',prob.Hosp_2_Death, ...
    'sigma',1/3, 'a',a, 'h',h, 'N',Polymod.UK_PP, 'beta',2.5, 'iota',0.6, 'N_total',total);

pi_val = 4.25e-6;

%% Build Age and Social Mixing Matrices
AgeMat0 = (Mix_H + Mix_S + Mix_W + Mix_O) .* (para.a .* para.h');
x = 10:-1:1;
epi = 0.3;
soc = SocioContMix(x, epi, total, D); 
soc = soc / sum(soc(:));

%% Contour Plots for Selected Pairs Including lambda_scale
paramNames = {'lambda_scale','rho', 'gamma', 'delta'};
nGrid = 50; scaling = 0.2;
pairs = nchoosek(1:length(paramNames), 2);
nPairs = size(pairs, 1);
nCols = 3; nRows = ceil(nPairs / nCols);

lambda_base = 1.0;

figure('Units','normalized','Position',[0.05 0.05 0.9 0.85]);
tiledlayout(nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');
labels = arrayfun(@(x) char('A'+x-1), 1:nPairs, 'UniformOutput', false);

% Greek label map
labels_map = containers.Map({'rho','gamma','delta','lambda_scale'}, ...
                             {'\rho','\gamma','\delta','\lambda'});

for k = 1:nPairs
    idx1 = pairs(k,1); idx2 = pairs(k,2);
    p1_name = paramNames{idx1}; p2_name = paramNames{idx2};

    % Set parameter ranges
    if strcmp(p1_name, 'lambda_scale')
        p1_vals = linspace(0.1, 1.0, nGrid);
    % elseif strcmp(p1_name, 'gamma')
    %     p1_vals = linspace(1/21, 1/7, nGrid);
    elseif strcmp(p1_name, 'delta')
        p1_vals = linspace(0.2, 1.0, nGrid);
    else
        p1_vals = linspace(para.(p1_name)*(1-scaling), para.(p1_name)*(1+scaling), nGrid);
    end

    if strcmp(p2_name, 'lambda_scale')
        p2_vals = linspace(0.1, 1.0, nGrid);
    % elseif strcmp(p2_name, 'gamma')
    %     p2_vals = linspace(1/21, 1/7, nGrid);
    elseif strcmp(p2_name, 'delta')
        p2_vals = linspace(0.2, 1.0, nGrid);
    else
        p2_vals = linspace(para.(p2_name)*(1-scaling), para.(p2_name)*(1+scaling), nGrid);
    end

    R0_grid = zeros(nGrid, nGrid);
    for i = 1:nGrid
        for j = 1:nGrid
            temp_para = para;
            lambda_eff = lambda_base;

            if strcmp(p1_name, 'lambda_scale')
                lambda_eff = p1_vals(i);
            elseif strcmp(p1_name, 'delta')
                temp_para.delta = para.delta * p1_vals(i);
            else
                temp_para.(p1_name) = p1_vals(i);
            end

            if strcmp(p2_name, 'lambda_scale')
                lambda_eff = p2_vals(j);
            elseif strcmp(p2_name, 'delta')
                temp_para.delta = para.delta * p2_vals(j);
            else
                temp_para.(p2_name) = p2_vals(j);
            end

            AgeMat_temp = AgeMat0 * lambda_eff;
            R0_grid(j,i) = computeR0_from_para(temp_para, pi_val, soc, AgeMat_temp);
        end
    end

    ax = nexttile;
    [X, Y] = meshgrid(p1_vals, p2_vals);
    contourf(X, Y, R0_grid, 20, 'LineColor', 'none');
    colormap(parula); colorbar;

    xlabel(['$', labels_map(p1_name), '$'], 'Interpreter', 'latex', 'FontSize', 20);
    ylabel(['$', labels_map(p2_name), '$'], 'Interpreter', 'latex', 'FontSize', 20);

    text(0.02, 0.95, labels{k}, 'Units','normalized', 'FontWeight','bold', 'FontSize',18, ...
        'BackgroundColor','w', 'Margin',1, 'EdgeColor','k', 'LineWidth',1.5, 'Parent',ax);
end

sgtitle('Contour Plots of $R_0$ vs Selected Parameter Pairs (using $\lambda$)', ...
    'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'latex');
print(gcf, 'R0_Contour_SelectedPairs_Lambda', '-dpng', '-r300');


%% Function to Compute R0
function R0 = computeR0_from_para(para, pi_val, soc, AgeMat0)
    term1 = ((1 - pi_val) * para.rho) / para.gamma;
    term2 = ((1 - pi_val) * (1 - para.rho)) ./ (para.delta + para.gamma);
    term3 = (para.iota * pi_val) ./ (para.delta + para.gamma);
    scalar_sum = term1 + term2 + term3;

    n_decile = size(soc, 1);
    n_age = size(AgeMat0, 1);
    AS_0 = kron(soc, AgeMat0);

    riskgroups = zeros(n_decile, n_decile);
    for i = 1:n_decile
        for j = 1:n_decile
            rows = ((i-1)*n_age + 1):(i*n_age);
            cols = ((j-1)*n_age + 1):(j*n_age);
            block = para.beta .* (AS_0(rows, cols) .* diag(scalar_sum));
            riskgroups(i,j) = max(real(eig(block)));
        end
    end

    R0 = max(eig(riskgroups));
end