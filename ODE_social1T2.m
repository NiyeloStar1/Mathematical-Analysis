%% The code for my PhD Project, this contains the age and social mixing
% Using a ODE15s for this analysis. 

%% Input 
% para - structure of paraeter
%ICs - Structure of initial conditions

%Mix_H - Age ixing in the house
%Mix_W - Age ixing in at work
%Mix_S - Age ixing in school
%Mix_O - Age ixing in other places
%beta - Total ixing that occur in all the locations
%N - Population size
% n- number of age groups
% soc - Soc demographic groups mixing matrix
% x - vector of mixing parttern


%% Output

%Classes - array of number of individual in each compartment

%% FUNCTION DEFINITION

function [Classes,R0,R02] = ODE_social1T2(para,ICs,maxtime,Mix_H,Mix_W,Mix_S,Mix_O,n,soc,x)
%Set tolerance
opts = odeset('RelTol',1e-4,'AbsTol',1e-3); %tolerance level%tolerance level
     tspan = [0:1:maxtime]; %time span


%Shape of the initial condition
 All = reshape([ICs.S1'; ICs.E1'; ICs.A1'; ICs.J1';ICs.I1'; ICs.H1'; ICs.R1'; ICs.D1'; ICs.Ir1';...
                 ICs.S2'; ICs.E2'; ICs.A2'; ICs.J2';ICs.I2';ICs.H2'; ICs.R2'; ICs.D2'; ICs.Ir2'; ...
                 ICs.S3'; ICs.E3'; ICs.A3'; ICs.J3';ICs.I3';ICs.H3'; ICs.R3'; ICs.D3'; ICs.Ir3'; ...
                 ICs.S4'; ICs.E4'; ICs.A4'; ICs.J4';ICs.I4';ICs.H4'; ICs.R4'; ICs.D4'; ICs.Ir4'; ...
                 ICs.S5'; ICs.E5'; ICs.A5'; ICs.J5';ICs.I5';ICs.H5'; ICs.R5'; ICs.D5'; ICs.Ir5'; ...
                 ICs.S6'; ICs.E6'; ICs.A6'; ICs.J6';ICs.I6';ICs.H6'; ICs.R6'; ICs.D6'; ICs.Ir6'; ...
                 ICs.S7'; ICs.E7'; ICs.A7'; ICs.J7';ICs.I7';ICs.H7'; ICs.R7'; ICs.D7'; ICs.Ir7'; ...
                 ICs.S8'; ICs.E8'; ICs.A8'; ICs.J8';ICs.I8';ICs.H8'; ICs.R8'; ICs.D8'; ICs.Ir8'; ...
                 ICs.S9'; ICs.E9'; ICs.A9'; ICs.J9';ICs.I9';ICs.H9'; ICs.R9'; ICs.D9'; ICs.Ir9'; ...
                 ICs.S10'; ICs.E10'; ICs.A10'; ICs.J10';ICs.I10';ICs.H10';ICs.R10'; ICs.D10'; ICs.Ir10'], []*n,1);

    [t , pop] = ode15s(@(t,y)diff_socialmodel(t,y,para),tspan,All,opts);
     

 Classes = struct('S1',pop(:,1:1:n),'E1',pop(:,n+1:1:2*n),'A1',pop(:,(2*n)+1:1:3*n),'J1',pop(:,(3*n)+1:1:4*n),'I1',pop(:,(4*n)+1:1:5*n),'H1',pop(:,(5*n)+1:1:6*n), ...
          'R1',pop(:,(6*n)+1:1:7*n),'D1',pop(:,(7*n)+1:1:8*n),'Ir1',pop(:,(8*n)+1:1:9*n),'S2',pop(:,(9*n)+1:1:10*n), ...
          'E2',pop(:,(10*n)+1:1:11*n),'A2',pop(:,(11*n)+1:1:12*n),'J2',pop(:,(12*n)+1:1:13*n),'I2',pop(:,(13*n)+1:1:14*n),'H2',pop(:,(14*n)+1:1:15*n), ...
          'R2',pop(:,(15*n)+1:1:16*n),'D2',pop(:,(16*n)+1:1:17*n),'Ir2',pop(:,(17*n)+1:1:18*n),'S3', pop(:,(18*n)+1:1:19*n),'E3', pop(:,(19*n)+1:1:20*n), ...
          'A3',pop(:,(20*n)+1:1:21*n),'J3',pop(:,(21*n)+1:1:22*n),'I3',pop(:,(22*n)+1:1:23*n),'H3',pop(:,(23*n)+1:1:24*n),'R3',pop(:,(24*n)+1:1:25*n), ...
          'D3',pop(:,(25*n)+1:1:26*n),'Ir3',pop(:,(26*n)+1:1:27*n),'S4',pop(:,(27*n)+1:1:28*n),'E4', pop(:,(28*n)+1:1:29*n),'A4', pop(:,(29*n)+1:1:30*n), ...
          'J4',pop(:,(30*n)+1:1:31*n),'I4',pop(:,(31*n)+1:1:32*n),'H4',pop(:,(32*n)+1:1:33*n),'R4',pop(:,(33*n)+1:1:34*n),'D4',pop(:,(34*n)+1:1:35*n), ...
          'Ir4',pop(:,(35*n)+1:1:36*n),'S5',pop(:,(36*n)+1:1:37*n),'E5',pop(:,(37*n)+1:1:38*n),'A5', pop(:,(38*n)+1:1:39*n),'J5',pop(:,(39*n)+1:1:40*n), ...
          'I5',pop(:,(40*n)+1:1:41*n),'H5',pop(:,(41*n)+1:1:42*n),'R5',pop(:,(42*n)+1:1:43*n),'D5',pop(:,(43*n)+1:1:44*n),'Ir5',pop(:,(44*n)+1:1:45*n), ...
          'S6',pop(:,(45*n)+1:1:46*n),'E6',pop(:,(46*n)+1:1:47*n),'A6',pop(:,(47*n)+1:1:48*n),'J6', pop(:,(48*n)+1:1:49*n),'I6', pop(:,(49*n)+1:1:50*n), ...
          'H6',pop(:,(50*n)+1:1:51*n),'R6',pop(:,(51*n)+1:1:52*n),'D6',pop(:,(52*n)+1:1:53*n),'Ir6',pop(:,(53*n)+1:1:54*n),'S7',pop(:,(54*n)+1:1:55*n), ...
          'E7',pop(:,(55*n)+1:1:56*n),'A7',pop(:,(56*n)+1:1:57*n),'J7',pop(:,(57*n)+1:1:58*n),'I7', pop(:,(58*n)+1:1:59*n),'H7', pop(:,(59*n)+1:1:60*n), ...
          'R7',pop(:,(60*n)+1:1:61*n),'D7',pop(:,(61*n)+1:1:62*n),'Ir7',pop(:,(62*n)+1:1:63*n),'S8',pop(:,(63*n)+1:1:64*n),'E8',pop(:,(64*n)+1:1:65*n), ...
          'A8',pop(:,(65*n)+1:1:66*n),'J8',pop(:,(66*n)+1:1:67*n),'I8',pop(:,(67*n)+1:1:68*n),'H8', pop(:,(68*n)+1:1:69*n),'R8', pop(:,(69*n)+1:1:70*n), ...
          'D8',pop(:,(70*n)+1:1:71*n),'Ir8',pop(:,(71*n)+1:1:72*n),'S9',pop(:,(72*n)+1:1:73*n),'E9',pop(:,(73*n)+1:1:74*n),'A9',pop(:,(74*n)+1:1:75*n), ...
          'J9',pop(:,(75*n)+1:1:76*n),'I9',pop(:,(76*n)+1:1:77*n),'H9',pop(:,(77*n)+1:1:78*n),'R9', pop(:,(78*n)+1:1:79*n),'D9', pop(:,(79*n)+1:1:80*n), ...
          'Ir9',pop(:,(80*n)+1:1:81*n),'S10',pop(:,(81*n)+1:1:82*n),'E10',pop(:,(82*n)+1:1:83*n),'A10',pop(:,(83*n)+1:1:84*n),'J10',pop(:,(84*n)+1:1:85*n), ...
          'I10',pop(:,(85*n)+1:1:86*n),'H10',pop(:,(86*n)+1:1:87*n),'R10',pop(:,(87*n)+1:1:88*n),'D10', pop(:,(88*n)+1:1:89*n),'Ir10', pop(:,(89*n)+1:1:90*n), ...
          't',t);

  %% Set up the population change

 function dpop = diff_socialmodel(t, pop, para)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Part 1: Logistic curve & R0 Calculation (t = 0 mixing)
    y_logistic = GeneralisedRichard(maxtime);
    % Compute pi for each age group (assumes same for all age classes)
    pi_val = zeros(length(y_logistic),1);
    for m = 1:length(y_logistic)
         pi_val(m) = y_logistic(m) / sum(sum(para.N)); % para.N_all should be the full population matrix
    end
    
    % Define the initial AgeMixing matrix at t=0 using all mixing locations
    AgeMat0 = (Mix_H + Mix_W + Mix_S + Mix_O) .* (para.a .* para.h');
    AS_0 = kron(soc, AgeMat0);
    nBlock = 21;   % number of age classes per social group
    nGroups = 10;  % number of social groups
    
    % Build AgeMSoc0 using loops (10x10 cell array)
    AgeMSoc0 = cell(nGroups, nGroups);
    for i = 1:nGroups
       for j = 1:nGroups
          rows = ( (i-1)*nBlock+1 : i*nBlock );
          cols = ( (j-1)*nBlock+1 : j*nBlock );
          AgeMSoc0{i,j} = AS_0(rows, cols);
       end
    end
    
    % Compute riskgroups Grand R0 as in your code
    riskgroups = zeros(nGroups, nGroups);
    
    
    % SCALAR TERM COMPUTATION
        term1 = ((1 - pi_val(1)) * para.rho) / para.gamma;  % Constant over age
        term2 = ((1 - pi_val(1)) * (1 - para.rho)) ./ (para.delta + para.gamma); % Vector of size 21
        term3 = (para.iota * pi_val(1)) ./ (para.delta + para.gamma);            % Vector of size 21

        scalar_sum = (term1 + term2 + term3); % Total scalar

    for i = 1:nGroups
       for j = 1:nGroups
          % Multiply each cell block by para.beta and the scalar transmission factor Sp.
        block_matrix = para.beta * AgeMSoc0{i,j} .* scalar_sum;
        % Compute its eigenvalues and take the maximum (using the real part)
        ev = eig(block_matrix);
        maxeig{i,j} = max(real(ev));
        riskgroups(i,j) = maxeig{i,j};  % Build the decile-level "risk groups" matrix.
       end
    end
    R0 = riskgroups ;
    R02 = max(eig(R0));
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Part 2: Determine current mixing based on time t

    ageMix = Mix_H + Mix_S + Mix_W + Mix_O;
         w = soc;
    % if (t > 0) && (t <= 83)
    %      ageMix = Mix_H + Mix_S + Mix_W + Mix_O;
    %      w = soc;
    % elseif (t > 83) && (t <= 164)
    %      ageMix = NPIchange(Mix_H, Mix_S, Mix_W, Mix_O, 'TOTAL_LOCKDOWN');
    %      w = 0.5 * soc;
    % elseif (t > 164) && (t <= 229)
    %      ageMix = NPIchange(Mix_H, Mix_S, Mix_W, Mix_O, 'EASING');
    %      w = soc;
    % elseif (t > 229) && (t <= 290)
    %      ageMix = NPIchange(Mix_H, Mix_S, Mix_W, Mix_O, 'RELAXED_RESTRICTION');
    %      w = soc;
    % else
    %      ageMix = NPIchange(Mix_H, Mix_S, Mix_W, Mix_O, 'RELAXED_RESTRICTION');
    %      w = soc;
    % end
    % 
    % Current mixing matrix
    AgeMat = ageMix .* (para.a .* para.h');
    AS = kron(w, AgeMat);
    % Build AgeSoc as a 10x10 cell array (each cell: nBlock x nBlock)
    AgeSoc = cell(nGroups, nGroups);
    for i = 1:nGroups
       for j = 1:nGroups
          rows = ((i-1)*nBlock+1 : i*nBlock);
          cols = ((j-1)*nBlock+1 : j*nBlock);
          AgeSoc{i,j} = AS(rows, cols);
       end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Part 3: Compute population derivatives for each social group
    % Total compartments per group = 10: S, E, A, J, Q, I, H, R, D, Ir
    totalCompartments = 9;
    n = 21; % number of age classes in each compartment
    dpop = zeros(length(pop),1);
    
    % Loop over each social group (i = 1,...,10)
    for i = 1:nGroups
       % Indices for group i in the pop vector
       idx_start = (i-1) * totalCompartments * n + 1;
       idx_end = i * totalCompartments * n;
       group_pop = pop(idx_start:idx_end);
       
       % Extract compartments for group i
       S = group_pop(1:n);
       E = group_pop(n+1:2*n);
       A = group_pop(2*n+1:3*n);
       J = group_pop(3*n+1:4*n);
       I = group_pop(4*n+1:5*n);
       H = group_pop(5*n+1:6*n);
       R = group_pop(6*n+1:7*n);
       D = group_pop(7*n+1:8*n);
       Ir = group_pop(8*n+1:9*n);
       % Ir = group_pop(9*n+1:10*n);
       
       % Compute SocInf for group i by summing contributions from all groups j
       SocInf = zeros(n,1);
       for j = 1:nGroups
          % For group j, extract the relevant compartments
          idx_j = (j-1) * totalCompartments * n + 1;
          group_pop_j = pop(idx_j : idx_j + totalCompartments*n - 1);
          A_j = group_pop_j(2*n+1:3*n);
          J_j = group_pop_j(3*n+1:4*n);
          I_j = group_pop_j(4*n+1:5*n);
          % I_j = group_pop_j(5*n+1:6*n);
          % Sum contributions using the precomputed AgeSoc cell
          SocInf = SocInf + AgeSoc{i,j} * (J_j +  para.iota *I_j +  A_j);
       end
       % Multiply by S and normalize by group population (para.N{i} should be a column vector of length n)
       SocInf = SocInf .* (S ./ para.D(i)');
       
       % Compute beta for group i
       beta_i = para.beta * SocInf;
       
       % Compute ODEs for group i (exactly as in your original code)
       dS = - beta_i + para.epsilon * R;
       dE = beta_i - para.sigma *((1 - pi_val(m)) * para.rho + (1 - pi_val(m)) * para.theta +pi_val(m))*E;
       dA = para.sigma*(1 - pi_val(m)) * para.rho * E - para.gamma * A;
       dJ = para.sigma*(1 - pi_val(m)) * para.theta * E - para.delta .* J - para.gamma * J;
       dI = para.sigma*pi_val(m) * E - para.delta .*I - para.gamma*I;
       dH = para.delta .* (J +  I) - para.d .* H - para.gamma * H;
       dR = para.gamma * (J +  A +  H  +  para.iota*I) - para.epsilon * R;
       dD = para.d .* H;
       dIr = para.sigma*pi_val(m) * E;
       % dIr = 1/para.PHI * para.p * Q;
       
       % Store the derivatives for group i
       group_dpop = [dS; dE; dA; dJ; dI; dH; dR; dD; dIr];
       dpop(idx_start:idx_end) = group_dpop;
    end
 end
end


