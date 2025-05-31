%% Run the ODE model
clc;close;clear

%% Set up the plot styling
% cmap = colormap(parula(10));
% set(gca, 'fontsize', 16)
% set(0, 'defaultaxeslinestyleorder', {'-*', ':', '-', '.'}) % or customize as needed
% set(0, 'defaultaxesfontsize', 16)
% set(0, 'defaultlinelinewidth', 2)

%% Load all required data
Polymod = load("MixingData.mat");
varyingpara = load("parameters.mat");
prob = load("Probab.mat");

% Load COVID-19 test data
Tdata = readtable('test_data.xlsx', 'Sheet', 'Jan_SepTest');

% Load the population data
Pop = load("EnglandDecilePop.mat");

N = Pop.TotalPop;
D = sum(N, 2)';  % Total by socioeconomic group
N_age = sum(N, 1);  % Total by age group
total = sum(D);  % Total England population

% Individual socio-demographic group data
N1 = N(1, :); N2 = N(2, :); N3 = N(3, :); N4 = N(4, :); N5 = N(5, :);
N6 = N(6, :); N7 = N(7, :); N8 = N(8, :); N9 = N(9, :); N10 = N(10, :);

%% Extract the age mixing matrices
Mix_H = Polymod.UK_from_toH;
Mix_S = Polymod.UK_from_toS;
Mix_W = Polymod.UK_from_toW;
Mix_O = Polymod.UK_from_toO;

%% Extract the infection parameter by age
a = varyingpara.Susceptibility; % Susceptibility
h = varyingpara.Detection; % Infectivity

%% Define the number of age and social groups
x = [10 9 8 7 6 5 4 3 2 1];  % Socio-demographic levels (1-10)
n = 21;  % Number of age groups
maxtime = 200;

%% Set up the parameters
para = struct('tau', 1, ...
    'epsilon', 2 / 100, ...
    'p', 0.8526, ...
    'rho', 0.6, 'theta', 0.4, ...
    'delta', prob.Sympt_2_hosp, ...
    'gamma', 1 / 14, 'gamma_2', 1 / 14, ...
    'd', prob.Hosp_2_Death, ...
    'sigma', 1 / 3, ...
    'a', a, 'h', h, 'N', Polymod.UK_PP, ...
    'beta', 2.5, ...
    'iota', 0.6, ...
    'N_total', total, 'N1', N1, 'N2', N2, 'N3', N3, 'N4', N4, ...
    'N5', N5, 'N6', N6, 'N7', N7, 'N8', N8, 'N9', N9, 'N10', N10, 'N_age', N_age, 'D', N, ...
    'N_social', D);

%% Define the initial conditions
EICs = load("ExposedIcs.mat");
E = EICs.EICs.Exposed;  % Extract exposed matrix
% I1 = [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
% I2 = E(2, :); I3 = E(3, :); I4 = E(4, :); I5 = E(5, :);
% I6 = E(6, :); I7 = E(7, :); I8 = E(8, :); I9 = E(9, :); I10 = E(10, :);


%% Vary epi and run ODE
epiVals = [0,0.3, 0.7, 1];  % Values for epi
nE = numel(epiVals);
ClassesAll = cell(nE, 1);

for k = 1:nE
    epi = epiVals(k);
    if epi == 0
        I1= [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
        I2 = I1; I3 = I1; I4 = I1; I5 = I1;
        I6 = I1; I7 = I1; I8 = I1; I9 = I1; I10 = I1 ;
    else
        I1 = [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
        I2 = E(2, :); I3 = E(3, :); I4 = E(4, :); I5 = E(5, :);
        I6 = E(6, :); I7 = E(7, :); I8 = E(8, :); I9 = E(9, :); I10 = E(10, :);
    end

    % Susceptible population
S1 = N1 - I1;

% Initial conditions struct
ICs = struct(...
    'S1', S1,    'E1', zeros(1,n), 'A1', zeros(1,n), 'J1', zeros(1,n), 'I1', I1,    'H1', zeros(1,n), 'R1', zeros(1,n), 'D1', zeros(1,n), 'Ir1', zeros(1,n), ...
    'S2', N2,    'E2', zeros(1,n), 'A2', zeros(1,n), 'J2', zeros(1,n), 'I2', I2,    'H2', zeros(1,n), 'R2', zeros(1,n), 'D2', zeros(1,n), 'Ir2', zeros(1,n), ...
    'S3', N3,    'E3', zeros(1,n), 'A3', zeros(1,n), 'J3', zeros(1,n), 'I3', I3,    'H3', zeros(1,n), 'R3', zeros(1,n), 'D3', zeros(1,n), 'Ir3', zeros(1,n), ...
    'S4', N4,    'E4', zeros(1,n), 'A4', zeros(1,n), 'J4', zeros(1,n), 'I4', I4,    'H4', zeros(1,n), 'R4', zeros(1,n), 'D4', zeros(1,n), 'Ir4', zeros(1,n), ...
    'S5', N5,    'E5', zeros(1,n), 'A5', zeros(1,n), 'J5', zeros(1,n), 'I5', I5,    'H5', zeros(1,n), 'R5', zeros(1,n), 'D5', zeros(1,n), 'Ir5', zeros(1,n), ...
    'S6', N6,    'E6', zeros(1,n), 'A6', zeros(1,n), 'J6', zeros(1,n), 'I6', I6,    'H6', zeros(1,n), 'R6', zeros(1,n), 'D6', zeros(1,n), 'Ir6', zeros(1,n), ...
    'S7', N7,    'E7', zeros(1,n), 'A7', zeros(1,n), 'J7', zeros(1,n), 'I7', I7,    'H7', zeros(1,n), 'R7', zeros(1,n), 'D7', zeros(1,n), 'Ir7', zeros(1,n), ...
    'S8', N8,    'E8', zeros(1,n), 'A8', zeros(1,n), 'J8', zeros(1,n), 'I8', I8,    'H8', zeros(1,n), 'R8', zeros(1,n), 'D8', zeros(1,n), 'Ir8', zeros(1,n), ...
    'S9', N9,    'E9', zeros(1,n), 'A9', zeros(1,n), 'J9', zeros(1,n), 'I9', I9,    'H9', zeros(1,n), 'R9', zeros(1,n), 'D9', zeros(1,n), 'Ir9', zeros(1,n), ...
    'S10', N10,  'E10', zeros(1,n), 'A10', zeros(1,n), 'J10', zeros(1,n), 'I10', I10, 'H10', zeros(1,n), 'R10', zeros(1,n), 'D10', zeros(1,n), 'Ir10', zeros(1,n) );

    soc = SocioContMix(x, epi, total, D);
    soc = soc / sum(soc(:));
    [Classes, ~, ~] = ODE_social1T2(para, ICs, maxtime, Mix_H, Mix_W, Mix_S, Mix_O, n, soc);
    ClassesAll{k} = Classes;
end

%% Plotting
colors = [
    0.64, 0.08, 0.18;   % Red
    0.00, 0.45, 0.74;   % Blue
    0.85, 0.33, 0.10;   % Orange
    % 0.47, 0.67, 0.19;   % Green
    % 0.49, 0.18, 0.56;   % Purple
    0.93, 0.69, 0.13;   % Yellow
    % 0.30, 0.75, 0.93;   % Cyan
];


figure;
t = linspace(0, maxtime, size(ClassesAll{1}.I1, 1));  % Time vector
lineStyles = {'-', '--', ':','-.'};  % Line styles for each epi value
% Prepare legend entries
legendEntries = arrayfun(@(x) sprintf('$\\epsilon = %.2f$', x), epiVals, 'UniformOutput', false);


for iGroup = 1:10
    subplot(2, 5, iGroup);
    hold on;
    for k = 1:nE
        Classes = ClassesAll{k};
        I_field = Classes.(['I' num2str(iGroup)]);
        I_tot = sum(I_field, 2);
        plot(t, I_tot, lineStyles{k}, 'LineWidth', 3.5,'Color', colors(k,:));
    end
    hold off;
    title(sprintf('Decile %d', iGroup));
    xlabel('\textbf{Time (days)}', 'Interpreter', 'latex',FontSize=18);
    ylabel('\textbf{Detected Infectious}', 'Interpreter', 'latex',FontSize=18);
    % if iGroup == 1
    % Add LaTeX legend to every subplot
    legend(legendEntries, 'Interpreter', 'latex', 'Location', 'northeast');
        % legend('$\\epsilon = %d$', (epiVals), 'Interpreter', 'latex', 'Location', 'best');
    % end
    grid on;
end

sgtitle('\textbf{Impact of Mixing Scenario on Infected by Deprivation Decile}', 'Interpreter', 'latex',FontSize=24);


