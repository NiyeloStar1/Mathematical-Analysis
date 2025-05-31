%% Run the ODE model
%Keep running tidy
clear all
close all
clc


%% Load all data needed at this point

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


epi = 1;
soc = SocioContMix(x,epi,total,D);

soc = soc/sum(sum(soc));
%% Define the number of age and social groups

% soc = 10;   %Sociodemographic group
n = 21;     %age group

%% Set up the Parameter

para = struct('epsilon',2/100, ...
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

maxtime = 200;

%% Define the initial conditions
EICs = load("ExposedIcs.mat"); %Load the exposed matrix
%Extract the exposed class
E = EICs.EICs.Exposed; %Exposed
I1 = [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
I2 = E(2,:);I3 =E(3,:);I4 =E(4,:);I5 =E(5,:);
I6 =E(6,:);I7 =E(7,:);I8 =E(8,:);I9 =E(9,:);I10 =E(10,:);

% Susceptible
S1 = N1 -I1;

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




%% Run the ODE
tic
[Classes,R0,R02] = ODE_social1T2(para, ICs, maxtime, Mix_H, Mix_W, Mix_S, Mix_O, n, soc);
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TOTAL PREVALENCE BY SOCIODEMOGRAPHIC GROUP %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GROUP ONE

prev1 = sum((Classes.A1 + Classes.J1 + Classes.I1 ),2);

%% GROUP TWO

prev2 = sum((Classes.A2 + Classes.J2 + Classes.I2 ),2);

%% GROUP THREE

prev3 = sum((Classes.A3 + Classes.J3 + Classes.I3 ),2);

%% GROUP FOUR

prev4 = sum((Classes.A4 + Classes.J4 +Classes.I4 ),2);

%% GROUP FIVE

prev5 = sum((Classes.A5 + Classes.J5 + Classes.I5 ),2);

%% GROUP SIX

prev6 = sum((Classes.A6 + Classes.J6  +Classes.I6 ),2);

%% GROUP SEVEN

prev7 = sum((Classes.A7 + Classes.J7  +Classes.I7 ),2);

%% GROUP EIGHT

prev8 = sum((Classes.A8 + Classes.J8 +Classes.I8 ),2);

%% GROUP NINE

prev9 = sum((Classes.A9 + Classes.J9 + Classes.I9 ),2);

%% GROUP TEN

prev10 = sum((Classes.A10 + Classes.J10 +Classes.I10 ),2);

PREV = [prev1,prev2,prev3,prev4,prev5,prev6,prev7,prev8,prev9,prev10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TOTAL INCIDENCE BY SOCIODEMOGRAPHIC GROUP %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IR = [sum(Classes.I1,2),sum(Classes.I2,2),sum(Classes.I3,2),sum(Classes.I4,2),...
    sum(Classes.I5,2),sum(Classes.I6,2),sum(Classes.I7,2),sum(Classes.I8,2),...
    sum(Classes.I9,2),sum(Classes.I10,2)];
% writematrix(IR,'infection.xlsx', 'Sheet', 'MixingTypes1')



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% HOSPITALISED BY SOCIODEMOGRAPHIC GROUP %%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hosp = [sum(Classes.H1,2),sum(Classes.H2,2),sum(Classes.H3,2),sum(Classes.H4,2),...
    sum(Classes.H5,2),sum(Classes.H6,2),sum(Classes.H7,2),sum(Classes.H8,2),...
    sum(Classes.H9,2),sum(Classes.H10,2)];

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%% Death BY SOCIODEMOGRAPHIC GROUP %%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Death = [sum(Classes.D1,2),sum(Classes.D2,2),sum(Classes.D3,2),sum(Classes.D4,2),...
    sum(Classes.D5,2),sum(Classes.D6,2),sum(Classes.D7,2),sum(Classes.D8,2),...
    sum(Classes.D9,2),sum(Classes.D10,2)];

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%DAILY  NUMBER OF DEATH REPORTED PER DAY BY SOCIODEMOGRAPHIC GROUP %%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
%% Group One

De1 = [zeros(1,n);diff(Classes.D1)];
Total_De1 = sum(De1,2);

%% Group Two

De2 = [zeros(1,n);diff(Classes.D2)];
Total_De2 = sum(De2,2);

%% Group Three

De3 = [zeros(1,n);diff(Classes.D3)];
Total_De3 = sum(De3,2);

%% Group Four

De4 = [zeros(1,n);diff(Classes.D4)];
Total_De4 = sum(De4,2);

%% Group Five

De5 = [zeros(1,n);diff(Classes.D5)];
Total_De5 = sum(De5,2);

%% Group Six

De6 = [zeros(1,n);diff(Classes.D6)];
Total_De6 = sum(De6,2);

%% Group Seven

De7 = [zeros(1,n);diff(Classes.D7)];
Total_De7 = sum(De7,2);

%% Group Eight

De8 = [zeros(1,n);diff(Classes.D8)];
Total_De8 = sum(De8,2);

%% Group Nine

De9 = [zeros(1,n);diff(Classes.D9)];
Total_De9 = sum(De9,2);

%% Group Ten

De10 = [zeros(1,n);diff(Classes.D10)];
Total_De10 = sum(De10,2);

TOTAL_DE = [Total_De1,Total_De2,Total_De3,Total_De4,Total_De5,...
    Total_De6,Total_De7,Total_De8,Total_De9,Total_De10];


%% Define the decile labels
decileLabels = { ...
    '1 – Most deprived', '2', '3', '4', ...
    '5', '6', '7', '8', '9', '10 – Least deprived'};

% Define 10 distinct colors and line styles
colors = lines(10);  % Generates 10 distinct colors
lineStyles = {'-', '--', ':', '-.', '-', '--', ':', '-.', '-', '--'};

xValues = 0:7;

figure;

% Subplot 1: Total Infection
subplot(2, 2, 1);
hold on;
handles = gobjects(10, 1);  % To store plot handles for legend
for i = 1:10
    handles(i) = plot(Classes.t / 30, PREV(:, i), ...
        'Color', colors(i, :), 'LineWidth', 3, 'LineStyle', lineStyles{i});
end
xticks(xValues);
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('Total Infection', 'Interpreter', 'latex', 'FontSize', 24);
xline(84 / 30, 'Color', 'r', 'LineWidth', 3);  % Lockdown Period
xline(164 / 30, 'Color', 'g', 'LineWidth', 3); % Lockdown Easing
legend(handles, decileLabels, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 18);
ax = gca;
ax.FontSize = 16;   % Increase tick label font size
hold off;

% Subplot 2: Confirmed Cases
subplot(2, 2, 2);
hold on;
for i = 1:10
    plot(Classes.t / 30, IR(:, i), ...
        'Color', colors(i, :), 'LineWidth', 3, 'LineStyle', lineStyles{i});
end
xticks(xValues);
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('Confirmed Cases', 'Interpreter', 'latex', 'FontSize', 24);
xline(84 / 30, 'Color', 'r', 'LineWidth', 3);
xline(164 / 30, 'Color', 'g', 'LineWidth', 3);
legend(handles, decileLabels, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 18);
ax = gca;
ax.FontSize = 16;
hold off;

% Subplot 3: Hospitalised
subplot(2, 2, 3);
hold on;
for i = 1:10
    plot(Classes.t / 30, Hosp(:, i), ...
        'Color', colors(i, :), 'LineWidth', 3, 'LineStyle', lineStyles{i});
end
xticks(xValues);
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('Hospitalised', 'Interpreter', 'latex', 'FontSize', 24);
xline(84 / 30, 'Color', 'r', 'LineWidth', 3);
xline(164 / 30, 'Color', 'g', 'LineWidth', 3);
legend(handles, decileLabels, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 12);
ax = gca;
ax.FontSize = 16;
hold off;

% Subplot 4: Daily Death
subplot(2, 2, 4);
hold on;
for i = 1:10
    plot(Classes.t / 30, TOTAL_DE(:, i), ...
        'Color', colors(i, :), 'LineWidth', 3, 'LineStyle', lineStyles{i});
end
xticks(xValues);
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('Daily Death', 'Interpreter', 'latex', 'FontSize', 24);
xline(84 / 30, 'Color', 'r', 'LineWidth', 3);
xline(164 / 30, 'Color', 'g', 'LineWidth', 3);
legend(handles, decileLabels, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 12);
ax = gca;
ax.FontSize = 16;
hold off;

% Add a shared title for the figure
sgtitle('\textbf{Disease Outcomes by Deprivation Decile. $\epsilon = 1$}', ...
    'Interpreter', 'latex', 'FontSize', 24);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%% HEAT MAPS %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Hospitalised
% 
Hospital = [sum(Classes.H1,1);sum(Classes.H2,1);sum(Classes.H3,1);sum(Classes.H4,1);...
    sum(Classes.H5,1);sum(Classes.H6,1);sum(Classes.H7,1);sum(Classes.H8,1);sum(Classes.H9,1);sum(Classes.H10,1)];



%% Infection
Infec = [sum(Classes.I1,1);sum(Classes.I2,1);sum(Classes.I3,1);sum(Classes.I4,1);...
    sum(Classes.I5,1);sum(Classes.I6,1);sum(Classes.I7,1);sum(Classes.I8,1);sum(Classes.I9,1);sum(Classes.I10,1)];

%% Death

CaseF = [sum(De1,1);sum(De2,1);sum(De3,1);sum(De4,1);sum(De5,1);sum(De6,1);sum(De7,1);sum(De8,1);sum(De9,1);sum(De10,1)];

%% Total Infection

%% GROUP ONE

prev1 = sum((Classes.A1 + Classes.J1 + Classes.I1 +Classes.H1),1);

%% GROUP TWO

prev2 = sum((Classes.A2 + Classes.J2 + Classes.I2 +Classes.H2),1);

%% GROUP THREE

prev3 = sum((Classes.A3 + Classes.J3 + Classes.I3 +Classes.H3),1);

%% GROUP FOUR

prev4 = sum((Classes.A4 + Classes.J4 + Classes.I4 +Classes.H4),1);

%% GROUP FIVE

prev5 = sum((Classes.A5 + Classes.J5 + Classes.I5 +Classes.H5),1);

%% GROUP SIX

prev6 = sum((Classes.A6 + Classes.J6 + Classes.I6 +Classes.H6),1);

%% GROUP SEVEN

prev7 = sum((Classes.A7 + Classes.J7 + Classes.I7 +Classes.H7),1);

%% GROUP EIGHT

prev8 = sum((Classes.A8 + Classes.J8 + Classes.I8 +Classes.H8),1);

%% GROUP NINE

prev9 = sum((Classes.A9 + Classes.J9 + Classes.I9 +Classes.H9),1);

%% GROUP TEN

prev10 = sum((Classes.A10 + Classes.J10 + Classes.I10 +Classes.H10),1);

PREV2 = [prev1;prev2;prev3;prev4;prev5;prev6;prev7;prev8;prev9;prev10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% HEAT MAP COMBINED PLOT %%%%%%%%%%%%%%%%%%%%%

customlabels2 = {'0-4', '5-9' ,'10-14', '15-19', '20-24', '25-29','30-34', ...
        '35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79', ...
        '80-84','85-89','90-94','95-99','99+'};
xValues = 1:21;
yValues = 1:10;

customlabels3 = {'1', '2' ,'3', '4', '5', '6','7','8','9','10'};

figure;

subplot(2,2,1)
clims1 = [min(PREV2(:)),max(PREV2(:))]; %Scale of the heatmap
imagesc(PREV2,clims1)
colormap(flipud(hot))
colorbar
title("(a) \textbf{Total infection}", 'Interpreter', 'latex', 'FontSize', 20)
xticks(xValues);
yticks(yValues);
xticklabels(customlabels2);
yticklabels(customlabels3);
xlabel('Age Groups', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('Deprivation Decile', 'Interpreter', 'latex', 'FontSize', 20)
ax = gca;            % get current axes
ax.FontSize = 16;    % set font size for tick labels

subplot(2,2,2)
clims2 = [min(Infec(:)),max(Infec(:))];
imagesc(Infec,clims2)
colormap(flipud(hot))
colorbar
title("(b) \textbf{Detected infection}", 'Interpreter', 'latex', 'FontSize', 20)
xticks(xValues);
yticks(yValues);
xticklabels(customlabels2);
yticklabels(customlabels3);
xlabel('Age Groups', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('Deprivation Decile', 'Interpreter', 'latex', 'FontSize', 20)
ax = gca;
ax.FontSize = 16;

subplot(2,2,3)
clims3 = [min(Hospital(:)),max(Hospital(:))];
imagesc(Hospital,clims3)
colormap(flipud(hot))
colorbar
title("(c) \textbf{Hospitalised}", 'Interpreter', 'latex', 'FontSize', 20)
xticks(xValues);
yticks(yValues);
xticklabels(customlabels2);
yticklabels(customlabels3);
xlabel('Age Groups', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('Deprivation Decile', 'Interpreter', 'latex', 'FontSize', 20)
ax = gca;
ax.FontSize = 16;

subplot(2,2,4)
clims4 = [min(CaseF(:)),max(CaseF(:))];
imagesc(CaseF,clims4)
colormap(flipud(hot))
colorbar
title("(d) \textbf{Death from Infection}", 'Interpreter', 'latex', 'FontSize', 20)
xticks(xValues);
yticks(yValues);
xticklabels(customlabels2);
yticklabels(customlabels3);
xlabel('Age Groups', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('Deprivation Decile', 'Interpreter', 'latex', 'FontSize', 20)
ax = gca;
ax.FontSize = 16;
