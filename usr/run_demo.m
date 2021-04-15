clear; close all;                % #ok<*NASGU> 

% set run parameters
runID    = 'demo3';              % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  20;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on (1) to display results
save_op  =  1;                   % switch on (1) to save output to file
plot_cv  =  1;                   % switch on (1) to live plot iterative convergence
bnchmrk  =  0;                   % switch on (1) to run manufactured solution benchmark

% set model domain parameters
L        =  40;                  % domain dimension
N        =  400 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  L/(N-2);             % grid spacing

% set model timing parameters
M        =  1e5;                 % number of time steps to take
tend     =  1e5;                 % end time for simulation [s]

% stress control parameters
Pu       = -1/20;                % ratio of pure-shear stress to buoyancy pressure
Si       =  0;                   % ratio of simple-shear stress to buoyancy pressure
B        =  1/2;                 % ratio of buoyancy pressure to tensile strength

% set model rheology parameters
n        =  0.75;                % non-Newtonian shear viscosity powerlaw
lambda   =  30;                  % exponential melt weakening factor
Es       =  10;                  % matrix viscosity activation energy
Ef       =  5;                   % melt viscosity activation energy
RHEAL    =  10;                  % damage healing rate
KDMG     =  10;                  % damage permeability enhancement factor
YDMG     =  0.2;                 % damage yield stress reduction factor
Pc       =  20;                  % deep crustal pressure (top boundary)

% thermo-chemical model parameters
Tc       =  0.20;                % deep crustal temperature (top boundary)
T0       =  0.50;                % initial temperature
T1       =  0.00;                % amplitude of random noise   
T2       =  0.01;                % amplitude of gaussian
D        =  0.4*L;               % thermal boundary layer thickness
PeT      =  1.0;                 % thermal Peclet number [u0.h0/kappaT0]
PeC      =  100;                 % major element Peclet number [u0.h0/Dc0]
St       =  3.0;                 % Stefan number [cp0.dT0/L0]
MAJ0     =  0.10;                % initial major element composition
MAJ1     =  0.05;                % amplitude of random noise   
MAJ2     = -0.0;                 % amplitude of gaussian
Da       =  100;                 % Dahmköhler number [t0/tr0]
PhDg     =  4.0;                 % Phase diagram scaling factor (> 1)
perCf    =  0.7;                 % peritectic liquidus composition
perCs    =  0.6;                 % peritectic solidus  composition
perT     =  0.4;                 % peritectic temperature
smx      =  (N/50)^2;            % smoothness of initial random noise in x-direction
smz      =  (N/50)^2;            % smoothness of initial random noise in z-direction
wx       =  L/5;                 % horizontal half-width of gaussian
wz       =  L/5;                 % vertical half-width of initial gaussian
xpos     =  0;                   % x-position of initial gaussian (0 = middle of domain)
zpos     = -L;                   % z-position of initial gaussian (0 = middle of domain)

% trace element model parameters
TRI0     =  1.0;                 % initial incompatible trace element composition
TRI1     =  0.0;                 % amplitude of random noise   
TRI2     =  0.0;                 % amplitude of gaussian
KTRI     =  0.01;                % incompatible trace element partitioning coefficient
TRC0     =  1.0;                 % initial compatible trace element composition
TRC1     =  0.0;                 % amplitude of random noise   
TRC2     =  0.0;                 % amplitude of gaussian
KTRC     =  100;                 % compatible trace element partitioning coefficient

% radiogenic isotope model parameters
IRP0     =  1.0;                 % initial radiogenic parent isotope composition
IRP1     =  0.0;                 % amplitude of random noise   
IRP2     =  1.0;                 % amplitude of gaussian
IRD0     =  0.01;                % initial radiogenic daughter isotope composition
IRD1     =  0.0;                 % amplitude of random noise   
IRD2     =  1.0;                 % amplitude of gaussian
DIRP     =  1.0;                 % parent isotope decay number (dimensionless half-life)
KIRP     =  10.;                 % parent isotope partitioning coefficient
KIRD     =  0.1;                 % daughter isotope partitioning coefficient

% stable isotope model parameters
ISS0     =  1.0;                 % initial solid stable isotope composition
ISS1     =  0.0;                 % amplitude of random noise   
ISS2     =  1.0;                 % amplitude of gaussian
ISF0     =  1.0;                 % initial fluid stable isotope composition
ISF1     =  0.0;                 % amplitude of random noise   
ISF2     =  1.0;                 % amplitude of gaussian

% set numerical model parameters
nup      =  50;                  % nonlinear coefficients, residual norms updated every 'nup' iterations
CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
theta    =  0.50;                % time-stepping scheme selector (1=BE, 1/2=CN, 0=FE)
rtol     =  1e-6;                % outer its relative tolerance
atol     =  1e-6;                % outer its absolute tolerance
minit    =  500;                 % minimum solver iterations
maxit    =  5e3;                 % maximum solver iterations
alpha    =  0.95;                % inner its step size (fraction of stable step) [0,1]
beta     =  0.80;                % iterative damping parameter (fraction of previous step) [0,1]
gamma    =  0.99;                % iterative relaxation for rheology updates [0,1]
kappa    =  1.00;                % regularisation of eIIvp for failure [0,1]
etamin   =  1e-2;                % minimum viscosity for regularisation
etamax   =  1e+3;                % maximum viscosity for stability
flim     =  1e-3;                % limit melt fraction in coefficients for stability

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
