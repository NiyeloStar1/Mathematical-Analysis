Polymod = load("MixingData.mat");
varyingpara = load("parameters.mat");
prob = load("Probab.mat");

% Load COVID-19 test data
Tdata = readtable('test_data.xlsx','Sheet','Jan_SepTest');

%% Load the population data 

Pop= load("EnglandDecilePop.mat");

N =  Pop.TotalPop;
D = (sum(N,2))';        %Total by socioeconomic group
N_age = (sum(N,1));     %Total by Age group
total = sum(D);         %Total England Population

N1= N(1,:); N2= N(2,:);N3= N(3,:);N4= N(4,:);N5= N(5,:);
N6= N(6,:); N7= N(7,:);N8= N(8,:);N9= N(9,:);N10= N(10,:);

%% Extract the age mixing matrix

Mix_H = Polymod.UK_from_toH;
Mix_S = Polymod.UK_from_toS;
Mix_W = Polymod.UK_from_toW;
Mix_O = Polymod.UK_from_toO;

%% Extract the infection parameter by Age

a = varyingpara.Susceptibility; % Susceptibility
h = varyingpara.Detection; % Infectivity

%% Setup the social mixing matrix
x = [10 9 8 7 6 5 4 3 2 1];


epi = 0.3;
soc = SocioContMix(x,epi,total,D);

soc = soc/sum(sum(soc));
%% Define the number of age and social groups

% soc = 10;   %Sociodemographic group
n = 21;     %age group

%% Set up the Parameter

para = struct('epsilon',2/100, ...
    'p',0.8526, ...
    'rho',0.6,'theta',0.4,...
    'delta',mean(prob.Sympt_2_hosp), ...
    'gamma',1/14,'gamma_2',1/14, ...
    'd',prob.Hosp_2_Death, ...
    'sigma',1/3,...
    'a',a,'h',h,'N',Polymod.UK_PP,...
    'beta',3.50,...%0.099
    'iota',0.6,...
    'N_total',total,'N1',N1,'N2',N2,'N3',N3,'N4',N4, ...
    'N5',N5,'N6',N6,'N7',N7,'N8',N8,'N9',N9,'N10',N10,'N_age',N_age,'D',N, ...
    'N_social',D);


%%%%%%%% CALCULATE R0
% ageMix = Mix_O+Mix_W+Mix_S+Mix_H;


% for i =1:10
% R0(i,:) = sum(para.beta*((sum(soc(i,:))*(para.a*ageMix.*para.h)))*1/21);
% end
% R0_total = sum(R0,2);
% 


%% Set up the run-time
maxtime = 250;

%% Define the initial conditions
EICs = load("ExposedIcs.mat"); %Load the exposed matrix
%Extract the exposed class
E = EICs.EICs.Exposed; %Exposed
I1 = ones(1,n); %[0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
I2 = E(2,:);I3 =E(3,:);I4 =E(4,:);I5 =E(5,:);
I6 =E(6,:);I7 =E(7,:);I8 =E(8,:);I9 =E(9,:);I10 =E(10,:);
% I1 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% J1 = [0 0 0 0 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0];

% Susceptible
S1 = N1 -I1;
% S2 = N2 -E2; S3 = N3 -E3; S4 = N4 -E4; S5 = N5 -E5; 
% S6 = N6 -E6;S7 = N7 -E7; S8 = N8 -E8; S9 = N9 -E9; S10 = N10 -E10; 



ICs = struct(...
    'S1', S1,    'E1', zeros(1,n), 'A1', zeros(1,n), 'J1', zeros(1,n), 'I1', I1,    'H1', zeros(1,n), 'R1', zeros(1,n), 'D1', zeros(1,n), 'Ir1', zeros(1,n), ...
    'S2', N2,    'E2', zeros(1,n), 'A2', zeros(1,n), 'J2', zeros(1,n), 'I2', zeros(1,n),    'H2', zeros(1,n), 'R2', zeros(1,n), 'D2', zeros(1,n), 'Ir2', zeros(1,n), ...
    'S3', N3,    'E3', zeros(1,n), 'A3', zeros(1,n), 'J3', zeros(1,n), 'I3', zeros(1,n),    'H3', zeros(1,n), 'R3', zeros(1,n), 'D3', zeros(1,n), 'Ir3', zeros(1,n), ...
    'S4', N4,    'E4', zeros(1,n), 'A4', zeros(1,n), 'J4', zeros(1,n), 'I4', zeros(1,n),    'H4', zeros(1,n), 'R4', zeros(1,n), 'D4', zeros(1,n), 'Ir4', zeros(1,n), ...
    'S5', N5,    'E5', zeros(1,n), 'A5', zeros(1,n), 'J5', zeros(1,n), 'I5', zeros(1,n),    'H5', zeros(1,n), 'R5', zeros(1,n), 'D5', zeros(1,n), 'Ir5', zeros(1,n), ...
    'S6', N6,    'E6', zeros(1,n), 'A6', zeros(1,n), 'J6', zeros(1,n), 'I6', zeros(1,n),    'H6', zeros(1,n), 'R6', zeros(1,n), 'D6', zeros(1,n), 'Ir6', zeros(1,n), ...
    'S7', N7,    'E7', zeros(1,n), 'A7', zeros(1,n), 'J7', zeros(1,n), 'I7', zeros(1,n),    'H7', zeros(1,n), 'R7', zeros(1,n), 'D7', zeros(1,n), 'Ir7', zeros(1,n), ...
    'S8', N8,    'E8', zeros(1,n), 'A8', zeros(1,n), 'J8', zeros(1,n), 'I8', zeros(1,n),    'H8', zeros(1,n), 'R8', zeros(1,n), 'D8', zeros(1,n), 'Ir8', zeros(1,n), ...
    'S9', N9,    'E9', zeros(1,n), 'A9', zeros(1,n), 'J9', zeros(1,n), 'I9', zeros(1,n),    'H9', zeros(1,n), 'R9', zeros(1,n), 'D9', zeros(1,n), 'Ir9', zeros(1,n), ...
    'S10', N10,  'E10', zeros(1,n), 'A10', zeros(1,n), 'J10', zeros(1,n), 'I10', zeros(1,n), 'H10', zeros(1,n), 'R10', zeros(1,n), 'D10', zeros(1,n), 'Ir10', zeros(1,n) );




% %% Run the ODE
tic
[Classes,R0,R02] = ODE_social1T2(para,ICs,maxtime,Mix_H,Mix_W,Mix_S,Mix_O,n,soc,x);
toc
% 
t = Classes.t;
% Total number of age groups
n = 21;

% Define custom vibrant colors (RGB triplets)
colors = [
    0.00, 0.45, 0.74;  % Vibrant blue     - S
    0.85, 0.33, 0.10;  % Deep orange      - E
    0.47, 0.67, 0.19;  % Bright green     - A
    0.49, 0.18, 0.56;  % Purple           - J
    0.93, 0.69, 0.13;  % Golden yellow    - I
    0.30, 0.75, 0.93;  % Cyan-blue        - H
    0.64, 0.08, 0.18;  % Crimson red      - R
];

compartments = {'S','E','A','J','I','H','R'};

%% Figure 1: Dynamics per Social Group (Summed Over Age Classes) — Publication Ready
% Uses your 'colors' and 'compartments' variables

fig1 = figure; clf;
set(fig1, 'Units','inches', 'Position',[0.6 0.6 14 10]);  % big canvas
t1 = tiledlayout(5, 2, 'Padding','compact', 'TileSpacing','compact');

lw_curve = 3.0;   % line thickness
lw_axes  = 1.6;   % axis spine thickness

axList1 = gobjects(10,1);
globalY1 = 0;

% Pre-create legend handles (from the first tile)
hLegend1 = gobjects(numel(compartments),1);

for g = 1:10
    % Sum over age groups for each compartment
    S = sum(Classes.(['S' num2str(g)]), 2);
    E = sum(Classes.(['E' num2str(g)]), 2);
    A = sum(Classes.(['A' num2str(g)]), 2);
    J = sum(Classes.(['J' num2str(g)]), 2);
    I = sum(Classes.(['I' num2str(g)]), 2);
    H = sum(Classes.(['H' num2str(g)]), 2);
    R = sum(Classes.(['R' num2str(g)]), 2);

    ax = nexttile(t1); axList1(g) = ax; hold(ax, 'on'); box(ax, 'on');
    set(ax, 'LineWidth', lw_axes, 'TickDir','out', ...
            'FontSize', 16, 'FontWeight','bold');

    % Plot lines (store handles in first tile for shared legend)
    if g == 1
        hLegend1(1) = plot(t, S, '-', 'LineWidth', lw_curve, 'Color', colors(1,:));
        hLegend1(2) = plot(t, E, '-', 'LineWidth', lw_curve, 'Color', colors(2,:));
        hLegend1(3) = plot(t, A, '-', 'LineWidth', lw_curve, 'Color', colors(3,:));
        hLegend1(4) = plot(t, J, '-', 'LineWidth', lw_curve, 'Color', colors(4,:));
        hLegend1(5) = plot(t, I, '-', 'LineWidth', lw_curve, 'Color', colors(5,:));
        hLegend1(6) = plot(t, H, '-', 'LineWidth', lw_curve, 'Color', colors(6,:));
        hLegend1(7) = plot(t, R, '-', 'LineWidth', lw_curve, 'Color', colors(7,:));
    else
        plot(t, S, '-', 'LineWidth', lw_curve, 'Color', colors(1,:));
        plot(t, E, '-', 'LineWidth', lw_curve, 'Color', colors(2,:));
        plot(t, A, '-', 'LineWidth', lw_curve, 'Color', colors(3,:));
        plot(t, J, '-', 'LineWidth', lw_curve, 'Color', colors(4,:));
        plot(t, I, '-', 'LineWidth', lw_curve, 'Color', colors(5,:));
        plot(t, H, '-', 'LineWidth', lw_curve, 'Color', colors(6,:));
        plot(t, R, '-', 'LineWidth', lw_curve, 'Color', colors(7,:));
    end

    title(ax, ['Decile ' num2str(g)], 'FontWeight','bold', 'FontSize', 18);
    xlabel(ax, 'Time(days)', 'Interpreter','latex', 'FontSize', 22, 'FontWeight','bold');
    ylabel(ax, 'Population',       'Interpreter','latex', 'FontSize', 22, 'FontWeight','bold');

    grid(ax, 'on'); ax.MinorGridLineStyle = '-'; grid(ax,'minor');
    set(ax, 'LooseInset', max(get(ax,'TightInset'), 0.02));

    % Track global ymax for unified limits
    localMax = max([S;E;A;J;I;H;R], [], 'all');
    if ~isempty(localMax), globalY1 = max(globalY1, localMax); end
end

% Unified Y-limits across all tiles (+5% headroom)
if globalY1 > 0
    yMax1 = ceil(globalY1 * 1.05);
    arrayfun(@(a) set(a, 'YLim', [0 yMax1]), axList1);
end

% Shared legend (boxed) for all tiles
% Shared legend across all tiles (create from an axes, not the layout)
lgd1 = legend(axList1(1), hLegend1, compartments, ...
              'Interpreter','latex', 'Orientation','horizontal', ...
              'NumColumns',4, 'Box','on');
set(lgd1, 'FontSize', 20, 'EdgeColor','k', 'LineWidth', 2.5);
lgd1.Layout.Tile = 'south';   % place below the tiled layout

% Main title
sgtitle(t1, 'Dynamics per Social Group (Summed Over Age Classes)', ...
        'Interpreter','latex', 'FontSize', 26, 'FontWeight','bold');

% Borderless, sharp export
exportgraphics(fig1, 'dynamics_per_social_group.png', ...
               'Resolution', 1000, 'BackgroundColor','white');
exportgraphics(fig1, 'dynamics_per_social_group.pdf', ...
               'ContentType','vector', 'BackgroundColor','white');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2: Infectiousness per Social Group (Summed Over Age Classes) — Publication Ready
compartments2 = {'A','J','I'};
cols2 = colors([3, 4, 7], :);   % use your A, J, I colors

fig2 = figure; clf;
set(fig2, 'Units','inches', 'Position',[0.6 0.6 14 10]);
t2 = tiledlayout(5, 2, 'Padding','compact', 'TileSpacing','compact');

lw_curve = 3.0; lw_axes = 1.6;
axList2 = gobjects(10,1);
globalY2 = 0;

hLegend2 = gobjects(numel(compartments2),1);

for g = 1:10
    A = sum(Classes.(['A' num2str(g)]), 2);
    J = sum(Classes.(['J' num2str(g)]), 2);
    I = sum(Classes.(['I' num2str(g)]), 2);

    ax = nexttile(t2); axList2(g) = ax; hold(ax, 'on'); box(ax, 'on');
    set(ax, 'LineWidth', lw_axes, 'TickDir','out', ...
            'FontSize', 16, 'FontWeight','bold');

    if g == 1
        hLegend2(1) = plot(t, A, '-', 'LineWidth', lw_curve, 'Color', cols2(1,:));
        hLegend2(2) = plot(t, J, '-', 'LineWidth', lw_curve, 'Color', cols2(2,:));
        hLegend2(3) = plot(t, I, '-', 'LineWidth', lw_curve, 'Color', cols2(3,:));
    else
        plot(t, A, '-', 'LineWidth', lw_curve, 'Color', cols2(1,:));
        plot(t, J, '-', 'LineWidth', lw_curve, 'Color', cols2(2,:));
        plot(t, I, '-', 'LineWidth', lw_curve, 'Color', cols2(3,:));
    end

    title(ax, ['Decile ' num2str(g)], 'FontWeight','bold', 'FontSize', 18);
    xlabel(ax, 'Time(days)', 'Interpreter','latex', 'FontSize', 22, 'FontWeight','bold');
    ylabel(ax, 'Population',       'Interpreter','latex', 'FontSize', 22, 'FontWeight','bold');

    grid(ax, 'on'); ax.MinorGridLineStyle = '-'; grid(ax,'minor');
    set(ax, 'LooseInset', max(get(ax,'TightInset'), 0.02));

    localMax = max([A;J;I], [], 'all');
    if ~isempty(localMax), globalY2 = max(globalY2, localMax); end
end

% Unified Y-limits (+5%)
if globalY2 > 0
    yMax2 = ceil(globalY2 * 1.05);
    arrayfun(@(a) set(a, 'YLim', [0 yMax2]), axList2);
end

% Shared, boxed legend
lgd2 = legend(axList2(1), hLegend2, compartments2, ...
              'Interpreter','latex', 'Orientation','horizontal', ...
              'NumColumns',3, 'Box','on');
set(lgd2, 'FontSize', 20, 'EdgeColor','k', 'LineWidth', 2.5);
lgd2.Layout.Tile = 'south';


% Main title
sgtitle(t2, 'Infectiousness per Social Group (Summed Over Age Classes)', ...
        'Interpreter','latex', 'FontSize', 26, 'FontWeight','bold');

% Borderless, sharp export
exportgraphics(fig2, 'infectiousness_per_social_group.png', ...
               'Resolution', 1000, 'BackgroundColor','white');
exportgraphics(fig2, 'infectiousness_per_social_group.pdf', ...
               'ContentType','vector', 'BackgroundColor','white');



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%% ALWAY UNCOMMENT THE OTHER SECTION IF YOU WANT A CONTOUR PLOT.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%% For this case, use the mean value of delta%%%%%%%%%%%%%
% %%%%%%%
% 
% % Define baseline parameters
% baseline = para;
% 
% % Define percent variation
% variation = 0.2;
% 
% % Time vector (adjust according to actual output from ODE function)
% t = 1:maxtime;
% 
% % Initialize storage for results
% I1_results = struct();
% 
% % Parameters to vary
% param_names = {'beta', 'gamma', 'rho', 'delta'};
% titles = {'\beta', '\gamma', '\rho', '\delta'};
% 
% % Start plotting
% figure;
% for i = 1:length(param_names)
%     pname = param_names{i};
% 
%     % Create three versions: -20%, baseline, +20%
%     vals = [1 - variation, 1, 1 + variation] * baseline.(pname);
% 
%     hold on
%     subplot(2, 2, i); % Create subplot grid 2x2
% 
%     for j = 1:3
%         % Modify parameter
%         para = baseline;
%         para.(pname) = vals(j);
% 
%         % Run simulation
%         [Classes, ~, ~] = ODE_social1T2(para, ICs, maxtime, Mix_H, Mix_W, Mix_S, Mix_O, n, soc, x);
% 
%         % Extract I1 (adjust this line depending on how your ODE returns it)
%         I1_total = sum(Classes.I1, 2);  % Sum across age groups
% 
%         % Store result
%         label = {'-20%', 'Baseline', '+20%'};
%         plot(t, I1_total, 'DisplayName', label{j}, 'LineWidth', 1.5); hold on;
%     end
% 
%     xlabel('Time (days)');
%     ylabel('I_1 (Infectious in group 1)');
%     title(['Variation in ', titles{i}]);
%     legend show;
%     grid on;
% end
% 
% sgtitle('Sensitivity of I_1 to Parameter Variation');
% 
