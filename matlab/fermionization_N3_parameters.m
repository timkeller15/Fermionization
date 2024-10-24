%% Simulation parameters for N = 3. Everything in dimensionless units as detailed in publication.
% System
N = 3; % Number of immersed particles
Ngrid = 512; % Position grid points (power of 2 for split-step algorithm)
NBEC = 1e4; % Number of BEC particles
gBEC = 1; % BEC intra-particle interaction strength
posmax = 7.5; % Position grid range
LTG = 12; % Impurity box potential range
wall = 1e8; % Value of the box potential edges
% Time-Evolution
dt = 1e-3; % Time step
cutoff = 1e-5; % Energy convergence criterion
steps = 5e4; % Number of time steps
samples = 5e3; % Number of data points
xi = LTG/N*(1 + 0.5*(mod(N,2) - 1)); % Peak position for seeding the initial pinned state wave function
% Interaction Loops
gMIXarr = 0:.01:5; % Range of the impurity-BEC interaction strength. We checked that the system energy is symmetric with respect to the sign of gMIX also numerically and therefore only considering positive gMIX is sufficient.
interactions = 0:.1:100; % Range of the intra-impurity interaction strength