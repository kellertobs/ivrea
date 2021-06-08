clear; close all;                % #ok<*NASGU> 

% set run parameters
runID    = 'magma3';              % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      = 20;                   % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on (1) to display results
save_op  =  1;                   % switch on (1) to save output to file
plot_cv  =  1;                   % switch on (1) to live plot iterative convergence
bnchmrk  =  0;                   % switch on (1) to run manufactured solution benchmark
diseq    =  0;                   % switch on (1) to use disequilibrium formulation, else local phase equilibrium imposed

% set model domain parameters
L        =  5;                   % domain dimension
N        =  200 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  L/(N-2);             % grid spacing

% set model timing parameters
M        =  1e5;                 % number of time steps to take
tend     =  1e5;                 % end time for simulation [s]

% stress control parameters
Pu       =  0;                   % ratio of pure-shear stress to buoyancy pressure
Si       =  0;                   % ratio of simple-shear stress to buoyancy pressure
Ty       =  1;                   % ratio of tensile strength to buoyancy pressure 

% set model rheology parameters
De       =  1e4;                 % visco-elastic Deborah number
n        =  0.0;                 % non-Newtonian shear viscosity powerlaw
lambda   =  30;                  % exponential melt weakening factor
Es       =  5;                   % matrix viscosity activation energy
Ef       =  0;                   % melt viscosity activation energy
EMAJ     =  1;                   % matrix viscosity dependence on major element composition
RHEAL    =  1;                   % damage healing rate
KDMG     =  1;                   % damage permeability enhancement factor
YDMG     =  0.05;                % damage yield stress reduction factor
DMG0     =  0.05;                % initial random damage amplitude
T0_eta   =  0.50;                % rheology reference temperature
f0_eta   =  0.0;                 % rheology reference melt fraction
e0_eta   =  10;                  % rheology reference strain rate
Pc       =  40;                  % crustal pressure (top boundary)

% thermo-chemical model parameters
Tc       =  0.00;                % crustal temperature (top boundary)
T0       =  0.001;               % initial temperature
T1       =  0.0;                 % amplitude of random noise (rel)
T2       =  0.70;                % amplitude of gaussian (abs)
D        =  0.60*L;              % thermal boundary layer thickness
PeT      =  2.0;                 % thermal Peclet number [u0.h0/kappaT0]
PeC      =  200;                 % major element Peclet number [u0.h0/Dc0]
St       =  2.4;                 % Stefan number [cp0.dT0/L0]
isotherm_topbot = 1;             % isothermal top/bot boundaries, else insulating
isotherm_sides  = 1;             % isothermal side    boundaries, else insulating
MAJ0     =  0.20;                % initial major element composition
MAJ1     =  0.10;                % amplitude of random noise   
MAJ2     = -0.00;                % amplitude of gaussian
Da       =  500;                 % Dahmköhler number [t0/tr0]
PhDg     =  5.0;                 % Phase diagram scaling factor (> 1)
perCf    =  0.75;                % peritectic liquidus composition
perCs    =  0.7;                 % peritectic solidus  composition
perT     =  0.3;                 % peritectic temperature
clap     =  7.e-4;               % Clapeyron slope for P-dependence of melting T
smx      =  (N/50)^2;            % smoothness of initial random noise in x-direction
smz      =  (N/50)^2;            % smoothness of initial random noise in z-direction
wx       =  L/3;                 % horizontal half-width of gaussian
wz       =  L/3;                 % vertical half-width of initial gaussian
xpos     =  0;                   % x-position of initial gaussian (0 = middle of domain)
zpos     = -L/2;                 % z-position of initial gaussian (0 = middle of domain)

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
DIRD     =  0.5;                 % daughter isotope decay number (dimensionless half-life)
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
CFL      =  1.0;                 % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'FROMM';             % advection scheme: UPW2, UPW3, FROMM, FLXDIV
theta    =  0.50;                % time-stepping scheme selector (1=BE, 1/2=CN, 0=FE)
rtol     =  1e-6;                % outer its relative tolerance
atol     =  1e-6;                % outer its absolute tolerance
minit    =  5;                   % minimum solver iterations
maxit    =  10;                  % maximum solver iterations
alpha    =  0.00; 
gamma    =  0.00;                % iterative relaxation for rheology updates [0,1]
delta    =  1e-4;                % numerical compressibility for P-diagonal stabilisation
kappa    =  0.00;                % regularisation of eIIvp for failure [0,1]
etamin   =  1e-5;                % minimum viscosity for regularisation
etamax   =  1e+2;                % maximum viscosity for stability
flim     =  1e-4;                % limit melt fraction in coefficients for stability

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
