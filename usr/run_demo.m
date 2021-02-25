clear; close all;                % #ok<*NASGU> 

% set run parameters
runID    = 'demo';               % run identifier
restart  =  12;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  10;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on (1) to display results
save_op  =  1;                   % switch on (1) to save output to file
plot_cv  =  1;                   % switch on (1) to live plot iterative convergence
bnchmrk  =  0;                   % switch on (1) to run manufactured solution benchmark

% set model domain parameters
L        =  50;                  % domain dimension
N        =  500 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  L/(N-2);             % grid spacing

% set model timing parameters
M        =  1e3;                 % number of time steps to take
tend     =  1e3;                 % end time for simulation [s]

% set model rheology parameters
n        =  1;                   % non-Newtonian shear viscosity powerlaw
Es       =  0;                   % matrix viscosity activation energy
Ef       =  0;                   % melt viscosity activation energy

% stress control parameters
Pu       =  0;                   % ratio of pure-shear stress to buoyancy pressure
Si       =  0;                   % ratio of simple-shear stress to buoyancy pressure
B        =  1;                   % ratio of buoyancy pressure to tensile strength

% thermo-chemical model parameters
Tc       =  0.20;                % deep crustal temperature (top boundary)
T0       =  0.40;                % initial temperature
T1       = -0.01;                % amplitude of random noise   
T2       = -0.00;                % amplitude of gaussian
D        =  L/3;                 % thermal boundary layer thickness
PeT      =  10;                  % thermal Peclet number [u0.h0/kappaT0]
PeC      =  100;                 % major element Peclet number [u0.h0/Dc0]
St       =  5;                   % Stefan number [cp0.dT0/L0]
MAJ0     =  0.20;                % initial major element composition
MAJ1     =  0.10;                % amplitude of random noise   
MAJ2     =  0.00;                % amplitude of gaussian
TRC0     =  1;                   % initial trace element composition
TRC1     =  0.1;                 % amplitude of random noise   
TRC2     =  0.2;                 % amplitude of gaussian
KTRC     =  0.01;                % trace element partitioning coefficient
SIS0     =  1;                   % initial stable isotope composition
SIS1     =  0.1;                 % amplitude of random noise   
SIS2     =  0.2;                 % amplitude of gaussian
Da       =  200;                 % Dahmköhler number [t0/tr0]
PhDg     =  4.0;                 % Phase diagram scaling factor (> 1)
perCf    =  0.8;                 % peritectic liquidus composition
perCs    =  0.7;                 % peritectic solidus  composition
perT     =  0.3;                 % peritectic temperature
smx      =  (N/50)^2;            % smoothness of initial random noise in x-direction
smz      =  (N/50)^2;            % smoothness of initial random noise in z-direction
wx       =  L/5;                 % horizontal half-width of gaussian
wz       =  L/5;                 % vertical half-width of initial gaussian
xpos     =  0;                   % x-position of initial gaussian (0 = middle of domain)
zpos     = -L/2;                 % z-position of initial gaussian (0 = middle of domain)

% set numerical model parameters
nup      =  25;                  % nonlinear coefficients, residual norms updated every 'nup' iterations
CFL      =  0.75;                % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'FRM';               % advection scheme ('UPW2', 'UPW3', or 'FRM')
theta    =  0.50;                % time-stepping scheme selector (1=BE, 1/2=CN, 0=FE)
rtol     =  1e-6;                % outer its relative tolerance
atol     =  1e-6;                % outer its absolute tolerance
minit    =  500;                 % minimum solver iterations
maxit    =  5e3;                 % maximum solver iterations
alpha    =  0.90;                % inner its step size (fraction of stable step) [0,1]
beta     =  0.75;                % iterative damping parameter (fraction of previous step) [0,1]
delta    =  0.99;               % iterative relaxation for rheology updates [0,1]
kappa    =  4.0;                % regularisation of eIIvp for failure [0,1]
etamin   =  0.03;                % minimum viscosity for regularisation
etamax   =  1e+2;                % maximum viscosity for stability

% create output directory
if ~isfolder(['../out/',runID])
    mkdir(['../out/',runID]);
end

% save input parameters and runtime options (unless restarting)
if restart == 0 
    parfile = ['../out/',runID,'/',runID,'_par'];
    save(parfile);
end

% run code
addpath ../src
main
