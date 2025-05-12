function y = GeneralisedRichard(maxtime)

r = 0.2245; %0.3898;
p= 0.9736; %0.8397;
alpha = 0.0839;%0.9677;
K= 20580000; %41757000;


% Define the function for the differential equation
dCdt = @(t, C) r * (C^p) * (1 - ((C/K)^alpha));

% Define the time span
tspan = [0:1:maxtime];

% Define the initial condition
C0 = 287;   % Initial value of C

% Solve the differential equation
[t, y] = ode45(dCdt, tspan, C0);

end
