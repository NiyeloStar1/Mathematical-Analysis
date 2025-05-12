% -------------------------------------------------------------------------
% Load Data
% -------------------------------------------------------------------------
Polymod = load("MixingData.mat");          % Age-based mixing data (home, school, work, other)
varyingpara = load("parameters.mat");      % Contains age-specific susceptibility and detection
prob = load("Probab.mat");                 % Transition probabilities (e.g., hospital to death)
Tdata = readtable('test_data.xlsx','Sheet','Jan_SepTest');  % Daily test data
Pop = load("EnglandDecilePop.mat");        % Population by socioeconomic decile and age
EICs = load("ExposedIcs.mat");             % Initial exposed cases

% -------------------------------------------------------------------------
% Define Population Structure
% -------------------------------------------------------------------------
N = Pop.TotalPop;              % Matrix of population by decile (rows) and age (columns)
D = sum(N,2)';                 % Total population by socioeconomic group (row sum)
N_age = sum(N,1);             % Total population by age group (column sum)
total = sum(D);               % Total population in England
n = 21;                       % Number of age groups
x = [10 9 8 7 6 5 4 3 2 1];   % Reversed decile ordering for social mixing

% -------------------------------------------------------------------------
% Extract Decile-Specific Populations
% -------------------------------------------------------------------------
for i = 1:10
    eval(sprintf('N%d = N(%d,:);', i, i));  % Each N1...N10 is a 1×21 vector for age distribution
end

% -------------------------------------------------------------------------
% Mixing Matrices
% -------------------------------------------------------------------------
Mix_H = Polymod.UK_from_toH;  % Home
Mix_S = Polymod.UK_from_toS;  % School
Mix_W = Polymod.UK_from_toW;  % Work
Mix_O = Polymod.UK_from_toO;  % Other

% -------------------------------------------------------------------------
% Age-Based Parameters
% -------------------------------------------------------------------------
a = varyingpara.Susceptibility;  % Age-specific susceptibility
h = varyingpara.Detection;       % Age-specific detection

% -------------------------------------------------------------------------
% Construct Social Mixing Matrix
% -------------------------------------------------------------------------
epi = 0.3;
soc = SocioContMix(x, epi, total, D);  % Custom function to create decile contact matrix
soc = soc / sum(soc(:));               % Normalize matrix

% -------------------------------------------------------------------------
% Model Parameters Structure
% -------------------------------------------------------------------------
para = struct('tau',1, ...
    'epsilon',2/100, ...
    'p',0.8526, ...
    'rho',0.6,'theta',0.4,...
    'delta',prob.Sympt_2_hosp, ...
    'gamma',1/14,'gamma_2',1/14, ...
    'd',prob.Hosp_2_Death, ...
    'sigma',1/3,...
    'a',a,'h',h,'N',Polymod.UK_PP,...
    'beta',2.5,...%0.099
    'iota',0.6,...
    'N_total',total,'N1',N1,'N2',N2,'N3',N3,'N4',N4, ...
    'N5',N5,'N6',N6,'N7',N7,'N8',N8,'N9',N9,'N10',N10,'N_age',N_age,'D',N, ...
    'N_social',D);
para.delta = mean(para.delta);

% -------------------------------------------------------------------------
% Initial Conditions Setup for Each Decile
% -------------------------------------------------------------------------
E = EICs.EICs.Exposed;
I1 = [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];  % Example seed for decile 1

for i = 2:10
    eval(sprintf('I%d = E(%d,:);', i, i));  % Extract initial exposed for deciles 2-10
end

for i = 1:10
    if i == 1
        Si = eval(sprintf('N%d - I%d;', i, i));  % Initial susceptible (population - infected)
    else
        Si = eval(sprintf('N%d;', i));           % Other deciles start fully susceptible
    end
    eval(sprintf([ ...
        'ICs.S%d = Si; ICs.E%d = zeros(1,n); ICs.A%d = zeros(1,n); ' ...
        'ICs.J%d = zeros(1,n); ICs.I%d = I%d; ICs.H%d = zeros(1,n); ' ...
        'ICs.R%d = zeros(1,n); ICs.D%d = zeros(1,n); ICs.Ir%d = zeros(1,n);'], ...
        i,i,i,i,i,i,i,i,i,i));
end

% -------------------------------------------------------------------------
% Simulation Duration
% -------------------------------------------------------------------------
maxtime = 350;  % Number of days to simulate

% -------------------------------------------------------------------------
% Custom Plot Colors (for visibility and aesthetics)
% -------------------------------------------------------------------------
colors = [
    0.00, 0.45, 0.74;   % Blue
    0.85, 0.33, 0.10;   % Orange
    0.47, 0.67, 0.19;   % Green
    0.49, 0.18, 0.56;   % Purple
    0.93, 0.69, 0.13;   % Yellow
    0.30, 0.75, 0.93;   % Cyan
    0.64, 0.08, 0.18;   % Red
];

% -------------------------------------------------------------------------
% Parameters to Vary and Their Labels
% -------------------------------------------------------------------------
param_names = {'beta', 'gamma', 'rho', 'delta'};
param_labels = {'q', '\gamma', '\rho', '\delta'};

% -------------------------------------------------------------------------
% Plotting Setup
% -------------------------------------------------------------------------
figure;
for i = 1:length(param_names)
    pname = param_names{i};    % Parameter name (e.g., 'beta')
    plabel = param_labels{i};  % LaTeX label (e.g., '\beta')
    base_val = para.(pname);   % Base value

    % Test values: ±20% variation + baseline
    test_vals = [0.6 0.8, 1, 1.2 1.4] * base_val;

    % Create subplot for current parameter
    subplot(2, 2, i);
    hold on;

    for j = 1:length(test_vals)
        % Create a temporary parameter set
        ptmp = para;
        ptmp.(pname) = test_vals(j);

        % Run ODE model
        [Classes, ~, ~] = ODE_social1T2(ptmp, ICs, maxtime, ...
                                        Mix_H, Mix_W, Mix_S, Mix_O, ...
                                        n, soc, x);

        % Total infectious individuals in decile 1
        I1_total = sum(Classes.I1, 2);
        t = Classes.t;

        % Plot trajectory
        lbl = sprintf('%s = %.2g', plabel, test_vals(j));
        plot(t, I1_total, 'LineWidth', 2.5, ...   % <== Thicker line
             'DisplayName', lbl, 'Color', colors(j,:));
    end

    % Axes labels and styling
    xlabel('Time (days)', 'Interpreter', 'latex',FontSize=18);
    ylabel(' (Dectected Infectious)', 'Interpreter', 'latex',FontSize=12);
    title(['Impact of varying $', plabel, '$'], 'Interpreter', 'latex', 'FontSize', 18);
    legend('Location', 'northeast');
    % grid on;

    
end

% Main title
sgtitle('Parameter Variations on Detected Infectious Individuals (Decile 1)');
