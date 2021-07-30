clear; close all;                % #ok<*NASGU> 

% set run parameters
runID    = 'demo5';               % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
plot_op  =  50;                  % switch on (1) to display results and print figures to file
save_op  =  0;                   % switch on (1) to save output to file
plot_cv  =  0;                   % switch on (1) to live plot iterative convergence
bnchmrk  =  0;                   % switch on (1) to run manufactured solution benchmark
diseq    =  0;                   % switch on (1) to use disequilibrium formulation, else local phase equilibrium imposed
inject   =  1;                   % switch on (1) injection of melt/heat/chemistry along base

% set model domain parameters
L        =  50;                  % domain dimension
N        =  200 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  L/(N-2);             % grid spacing

% set model timing parameters
M        =  1e5;                 % number of time steps to take
tend     =  1e5;                 % end time for simulation [s]
dt       =  1e-5;

% stress control parameters
Pu       =  -0.1/L;              % ratio of pure-shear stress to buoyancy pressure 
Si       =  0;                   % ratio of simple-shear stress to buoyancy pressure
Ty       =  2;                   % ratio of tensile strength to buoyancy pressure 
 
% set model rheology parameters
De       =  6e3;                 % visco-elastic Deborah number
n        =  0.5;                 % non-Newtonian shear viscosity powerlaw
lambda   =  25;                  % exponential melt weakening factor
Es       =  4;                   % matrix viscosity activation energy
Ef       =  0;                   % melt viscosity activation energy
EMAJ     =  1e-7;                % matrix viscosity dependence on major element composition
RHEAL    =  2;                   % damage healing rate
KDMG     =  1;                   % damage permeability enhancement factor
YDMG     =  0.05;                % damage yield stress reduction factor
DMG0     =  0.05;                % initial random damage amplitude
T0_eta   =  0.55;                % rheology reference temperature
f0_eta   =  0.0;                 % rheology reference melt fraction
e0_eta   =  1;                   % rheology reference strain rate
MAJ0_eta =  0.10;                % rheology reference composition
Pc       =  50;                  % crustal pressure (top boundary)

% thermo-chemical model parameters
f0       =  0.02;                % reference melt fraction
T0       =  0.54;                % initial temperature
T1       =  0.0;                 % amplitude of random noise (rel)
T2       =  0.000;               % amplitude of gaussian (abs)
Tc       =  0.00;                % crustal temperature (top boundary)
DLAB     =  0.60*L;              % thermal boundary (LAB) depth
DMOHO    =  0.20*L;              % crust-mantle boundary (MOHO) depth
PeT      =  0.25;                % thermal Peclet number [u0.h0/kappaT0]
PeC      =  500;                 % major element Peclet number [u0.h0/Dc0]
St       =  2.5;                 % Stefan number [cp0.dT0/L0]
isotherm_topbot = 1;             % isothermal top/bot boundaries, else insulating
isotherm_sides  = 0;             % isothermal side    boundaries, else insulating
MAJ0     =  0.15;                % initial major element composition
MAJ1     =  0.05;                % amplitude of random noise   
MAJ2     = -0.05;                % amplitude of gaussian
MAJc     =  0.15;                % crustal composition (MOHO to top)
Da       =  1e-16;               % Dahmköhler number [t0/tr0]
PhDg     =  5.0;                 % Phase diagram scaling factor (> 1)
perCf    =  0.75;                % peritectic liquidus composition
perCs    =  0.65;                % peritectic solidus  composition
perT     =  0.3;                 % peritectic temperature
clap     =  5e-4;                % Clapeyron slope for P-dependence of melting T
dTH2O    =  2e-3;                % solidus shift from water content
smx      =  (N/50)^2;            % smoothness of initial random noise in x-direction
smz      =  (N/50)^2;            % smoothness of initial random noise in z-direction
wx       =  4*L;                 % horizontal half-width of gaussian
wz       =  L/4;                 % vertical half-width of initial gaussian
xpos     =  0;                   % x-position of initial gaussian (0 = middle of domain)
zpos     = -L/1;                 % z-position of initial gaussian (0 = middle of domain)

% trace element model parameters
TRI0     =  1.0;                 % initial incompatible trace element composition
TRI1     =  0.05;                % amplitude of random noise   
TRI2     = -0.30;                % amplitude of gaussian
KTRI     =  0.01;                % incompatible trace element partitioning coefficient
TRC0     =  1.0;                 % initial compatible trace element composition
TRC1     = -0.05;                % amplitude of random noise   
TRC2     =  0.10;                % amplitude of gaussian
KTRC     =  100;                 % compatible trace element partitioning coefficient

% radiogenic isotope model parameters
IRP0     =  1.0;                 % initial radiogenic parent isotope composition
IRP1     =  0.05;                % amplitude of random noise   
IRP2     =  1.0;                 % amplitude of gaussian
IRD0     =  0.01;                % initial radiogenic daughter isotope composition
IRD1     =  0.05;                % amplitude of random noise   
IRD2     =  1.0;                 % amplitude of gaussian
DIRP     =  2.0;                 % parent isotope decay number (dimensionless half-life)
DIRD     =  1.0;                 % daughter isotope decay number (dimensionless half-life)
KIRP     =  10.;                 % parent isotope partitioning coefficient
KIRD     =  0.1;                 % daughter isotope partitioning coefficient

% stable isotope model parameters
ISS0     =  1.0;                 % initial solid stable isotope composition
ISS1     =  0.05;                % amplitude of random noise   
ISS2     =  1.0;                 % amplitude of gaussian
ISF0     =  1.0;                 % initial fluid stable isotope composition
ISF1     =  0.05;                % amplitude of random noise   
ISF2     =  1.0;                 % amplitude of gaussian

% set numerical model parameters
CFL      =  0.90;                % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'FROMM';             % advection scheme: UPW2, UPW3, FROMM, FLXDIV
theta    =  0.50;                % time-stepping scheme selector (1=BE, 1/2=CN, 0=FE)
rtol     =  1e-3;                % outer its relative tolerance
atol     =  1e-7;                % outer its absolute tolerance
minit    =  1;                   % minimum solver iterations
maxit    =  10;                  % maximum solver iterations
alpha    =  0.60;                % iterative relaxation for thermochem compaction updates [0,1]
gamma    =  0.60;                % iterative relaxation for rheology updates [0,1]
delta    =  3e-4;                % numerical compressibility for P-diagonal stabilisation
kappa    =  0.00;                % regularisation of eIIvp for failure [0,1]
etamin   =  1e-2;                % minimum viscosity for regularisation
etamax   =  1e+4;                % maximum viscosity for stability
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
