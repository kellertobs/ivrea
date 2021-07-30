% print run header
fprintf('\n\n*****  ivrea  |  %s  |  %s  *****\n\n',runID,datetime);

% load custom colormap
load ocean;  

ic  =  2:N-1 ;
ip  =  2:N   ;
im  =  1:N-1 ;
ibz = [2 N-1];  % boundary indices for closed top/bot boundaries
if abs(Si)>abs(Pu)
    ibx = [N-1 2];  % boundary indices for periodic side boundaries
else
    ibx = [2 N-1];  % boundary indices for closed side boundaries
end

% produce smooth random perturbations
rng(5);
dr = randn(N,N);
for i = 1:max(smx,smz)
    dr(ic,ic) = dr(ic,ic) ...
              + smz./max(smx,smz).*diff(dr(:,ic),2,1)./8 ...
              + smx./max(smx,smz).*diff(dr(ic,:),2,2)./8;
    dr = dr - mean(dr(:)); 
    dr = dr./max(abs(dr(:)));
    dr([1 2 end-1 end],:) = 0;
    dr(:,[1 2 end-1 end]) = 0;
end

% get mapping arrays
NP =  N   * N  ;
NW = (N-1)* N  ;
NU =  N  *(N-1);
MapP = reshape(1:NP,N  ,N  );
MapW = reshape(1:NW,N-1,N  );
MapU = reshape(1:NU,N  ,N-1) + NW;

% get coordinate arrays
z     = -0  -h/2:h:L  +h/2;
x     = -L/2-h/2:h:L/2+h/2;
xc    = (x(im)+x(ip))./2;
zc    = (z(im)+z(ip))./2;
[X,Z] = meshgrid(x,z);

% set gaussian perturbation shape function
gs   = exp(-(X+xpos).^2./wx^2).*exp(-(Z+zpos).^2./wz^2);

% initialise background velocity fields
Pu  = Pu+1e-16;  Si = Si+1e-16;
WP  =    (Z(im,:)+Z(ip,:))/2 *Pu;  % pure   shear WBG
WS  =   -(X(im,:)+X(ip,:))/2 *Si;  % simple shear WBG
UP  =   -(X(:,im)+X(:,ip))/2 *Pu;  % pure   shear UBG
US  = (L-(Z(:,im)+Z(:,ip))/2)*Si;  % simple shear UBG
WS  = 0.*WS; US = 2.*US; 
WBG = WP+WS;
UBG = UP+US;     
clear  WP WS UP US

% initialise geochemical signature
TRI  =  TRI0 .* (1 + TRI1.*dr + TRI2.*gs); TRIin = TRI;  TRIf = TRI/KTRI; % incompatible trace
TRC  =  TRC0 .* (1 + TRC1.*dr + TRC2.*gs); TRCin = TRC;  % compatible trace
IRP  =  IRP0 .* (1 + IRP1.*dr + IRP2.*gs); IRPin = IRP;  % radiogenic parent isotope
IRD  =  IRP.*DIRD./DIRP;                   IRDin = IRD;  % daughter isotope in decay equilibrium
ISS  =  ISS0 .* (1 + ISS1.*dr + ISS2.*gs); ISSin = ISS;  % solid stable isotope
ISF  =  ISF0 .* (1 + ISF1.*dr + ISF2.*gs); ISFin = ISF;  % fluid stable isotope

% initialise thermo-chemical solution and residual fields
if bnchmrk
    mms;
    f   = f_mms;
    T   = T_mms;
    MAJ = C_mms;
else
    T   =  Tc-(T0-Tc).*erf((-Z)./DLAB).*(1 + T1.*dr) + T2.*gs; Tin = T;       % temperature
    MAJ =  MAJ0 .* (1 + MAJ1.*dr) + MAJ2.*gs; MAJin = MAJ;                 % major elements
    MAJ = MAJ + (MAJc-MAJ0) .* (1-erf(75*(Z-DMOHO)/L))/2;
    f   =  0.*T;  fq = f+1;
    Pt  =  Pc + 5*Z;                                                       % lithostatic pressure
    while norm(fq-f,2)/N > 1e-12
        [fq,MAJsq,MAJfq] = equilibrium(T,MAJ,Pt,TRIf,perT,perCs,perCf,clap,dTH2O,PhDg); % melt fraction at equilibrium
        f    = 0.9*f    + 0.1*fq;
        TRIf = 0.9*TRIf + 0.1*(TRI./(f + (1-f).*KTRI));
    end
    fin = fq;
end
fprintf('\n\n*****  initial condition: T0 = %1.3f;  C0 = %1.3f;  f0 = %1.3f;\n\n',T0,MAJ0,f0);
    
% update phase major element composition
KMAJ = MAJsq./MAJfq;
MAJf = MAJ./(f + (1-f).*KMAJ);
MAJs = MAJ./(f./KMAJ + (1-f));
MAJf(f<=  1e-3) = MAJfq(f<=1e-3);  MAJs(f<=  1e-3) = MAJ(f<=    1e-3);
MAJf(f>=1-1e-3) = MAJ(f>=1-1e-3);  MAJs(f>=1-1e-3) = MAJsq(f>=1-1e-3);
res_f   = 0.*f;
res_T   = 0.*T;
res_MAJ = 0.*MAJ;

% initialise mechanical solution and residual fields
W  =  1.*WBG;  FW = 0.*W;
U  =  1.*UBG;  FU = 0.*U;
P  =  0.*f;    FP = 0.*P;
u  =  0.*U;    Fu = u;
w  =  0.*W;    Fw = w;
p  =  0.*P;    Fp = p;  po = p;   Pe = Pt;
S  = [W(:);U(:);w(:);u(:);P(:);p(:)];
wf = W+w./max(flim,(f(im,:)+f(ip,:))./2);
uf = U+u./max(flim,(f(:,im)+f(:,ip))./2);
DMG  = DMG0*(1+dr);  res_DMG = 0.*DMG;

% initialise parameter fields
RctR_f =  0.*f;
Rf     =  0.*f;
RT     =  0.*T;
RMAJ   =  0.*MAJ;
RTRI   =  0.*TRI;
RTRC   =  0.*TRC;
RIRP   =  0.*IRP;
RIRD   =  0.*IRD;
RISS   =  0.*ISS;
RISF   =  0.*ISF;
RDMG   =  0.*DMG;
ISR    =  ISS;
ups    =  0.*P;  upss   =  0.*P;  upsf = 0.*P;  Div_fV = 0.*P(ic,ic);
eps0   =  max(1/L,abs(Pu) + abs(Si)) + 1e-16;  
exx    =  0.*P - Pu;  ezz = 0.*P + Pu;  exz = zeros(N-1,N-1) - Si;  eps = 0.*P + (abs(Pu) + abs(Si));
epsVIS =  0.*eps;  epsELA = eps; epsDMG = 0.*eps;
upsVIS =  0.*ups;  upsELA = ups; upsDMG = 0.*ups;
txx    =  0.*exx;  tzz = 0.*ezz;  txz = 0.*exz;  tau = 0.*eps;  tauo = tau;
chi_vep =  0.*eps;  xi_vep  = 0.*ups;

% initialise timing parameters
it     =  0;
step   =  0;
time   =  0;

% print initial condition
if ~restart; up2date; output; end

% overwrite fields from file if restarting run
if     restart < 0  % restart from last continuation frame
    if isfile(['../out/',runID,'/',runID,'_cont.mat'])
        name = ['../out/',runID,'/',runID,'_par'];
        load(name);
        name = ['../out/',runID,'/',runID,'_cont.mat'];
        load(name);
        % increment time step
        time = time+dt;
        step = step+1;
    else
        restart = 0;
    end
elseif restart > 0  % restart from specified continuation frame
    name = ['../out/',runID,'/',runID,'_',num2str(restart),'.mat'];
    load(name);
    name = ['../out/',runID,'/',runID,'_par'];
    load(name);
    % increment time step
    time = time+dt;
    step = step+1;
end

% time stepping loop
while time < tend && step < M
    
    fprintf(1,'\n\n*****  step %d;  dt = %4.4e;  time = %4.4e;\n\n',step,dt,time);
    tic;
    
    % store previous solution
    fo   = f;    Rfo   = Rf;
    To   = T;    RTo   = RT;
    MAJo = MAJ;  RMAJo = RMAJ;
    TRIo = TRI;  RTRIo = RTRI;
    TRCo = TRC;  RTRCo = RTRC;
    IRPo = IRP;  RIRPo = RIRP;
    IRDo = IRD;  RIRDo = RIRD;
    ISSo = ISS;  RISSo = RISS;
    ISFo = ISF;  RISFo = RISF;
    DMGo = DMG;  RDMGo = RDMG;
    Pto  = Pt;   TRIfo = TRIf;
    Div_fVo = Div_fV;
    txxo = txx; tzzo = tzz; txzo = txz;
    tauo = tau; % - flxdiv(tau,U,W,h,ADVN,'adv').*dt;
    po   = p;   % - flxdiv(p  ,U,W,h,ADVN,'adv').*dt;
    dto  = dt;
    
    % reset residual norms and iteration count
    resnorm  = 1e3;
    resnorm0 = resnorm;
    it       = 0;
    
    % non-linear iteration loop
    startup = 2*double(step<=1) + double(step>1);
    while resnorm/resnorm0 >= rtol/startup^2 && resnorm >= atol/startup && it < maxit*startup || it < minit
        
        % increment iteration count
        it = it+1;
        
        % update thermo-chemical evolution
        thermochem; 
        
        % update nonlinear coefficients & auxiliary fields
        up2date;
        
        % solve fluid-mechanics equations
        fluidmech;
        
    end
            
    % print diagnostics
    fprintf(1,'\n         time to solution = %4.4f sec\n\n',toc);

    fprintf(1,'         min U    = %s%4.4f;   mean U    = %s%4.4f;   max U    = %s%4.4f;\n'  ,int8(min(  U(:))<0),min(  U(:)),int8(mean(  U(:))<0),mean(  U(:)),int8(max(  U(:))<0),max(  U(:)));
    fprintf(1,'         min W    = %s%4.4f;   mean W    = %s%4.4f;   max W    = %s%4.4f;\n'  ,int8(min( -W(:))<0),min( -W(:)),int8(mean( -W(:))<0),mean( -W(:)),int8(max( -W(:))<0),max( -W(:)));
    fprintf(1,'         min P    = %s%4.4f;   mean P    = %s%4.4f;   max P    = %s%4.4f;\n\n',int8(min(  P(:))<0),min(  P(:)),int8(mean(  P(:))<0),mean(  P(:)),int8(max(  P(:))<0),max(  P(:)));

    fprintf(1,'         min u    = %s%4.4f;   mean u    = %s%4.4f;   max u    = %s%4.4f;\n'  ,int8(min(  u(:))<0),min(  u(:)),int8(mean(  u(:))<0),mean(  u(:)),int8(max(  f(:))<0),max(  u(:)));
    fprintf(1,'         min w    = %s%4.4f;   mean w    = %s%4.4f;   max w    = %s%4.4f;\n'  ,int8(min( -w(:))<0),min( -w(:)),int8(mean( -w(:))<0),mean( -w(:)),int8(max(  f(:))<0),max( -w(:)));
    fprintf(1,'         min p    = %s%4.4f;   mean p    = %s%4.4f;   max p    = %s%4.4f;\n\n',int8(min(  p(:))<0),min(  p(:)),int8(mean(  p(:))<0),mean(  p(:)),int8(max(  f(:))<0),max(  p(:)));
    
    fprintf(1,'         min f    = %s%4.4f;   mean f    = %s%4.4f;   max f    = %s%4.4f;\n'  ,int8(min(  f(:))<0),min(  f(:)),int8(mean(  f(:))<0),mean(  f(:)),int8(max(  f(:))<0),max(  f(:)));
    fprintf(1,'         min T    = %s%4.4f;   mean T    = %s%4.4f;   max T    = %s%4.4f;\n'  ,int8(min(  T(:))<0),min(  T(:)),int8(mean(  T(:))<0),mean(  T(:)),int8(max(  T(:))<0),max(  T(:)));
    fprintf(1,'         min C    = %s%4.4f;   mean C    = %s%4.4f;   max C    = %s%4.4f;\n\n',int8(min(MAJ(:))<0),min(MAJ(:)),int8(mean(MAJ(:))<0),mean(MAJ(:)),int8(max(MAJ(:))<0),max(MAJ(:)));

    fprintf(1,'         min K    = %s%4.4f;   mean K    = %s%4.4f;   max K    = %s%4.4f;\n'  ,int8(min(log10(   K(:)))<0),min(log10(   K(:))),int8(mean(log10(   K(:)))<0),mean(log10(   K(:))),int8(max(log10(   K(:)))<0),max(log10(   K(:))));
    fprintf(1,'         min  eta = %s%4.4f;   mean  eta = %s%4.4f;   max  eta = %s%4.4f;\n'  ,int8(min(log10( eta(:)))<0),min(log10( eta(:))),int8(mean(log10( eta(:)))<0),mean(log10( eta(:))),int8(max(log10( eta(:)))<0),max(log10( eta(:))));
    fprintf(1,'         min zeta = %s%4.4f;   mean zeta = %s%4.4f;   max zeta = %s%4.4f;\n\n',int8(min(log10(zeta(:)))<0),min(log10(zeta(:))),int8(mean(log10(zeta(:)))<0),mean(log10(zeta(:))),int8(max(log10(zeta(:)))<0),max(log10(zeta(:))));
    
    fprintf(1,'         min ups  = %s%4.4f;   mean ups  = %s%4.4f;   max ups  = %s%4.4f;\n'  ,int8(min(ups(:))<0),min(ups(:)),int8(mean(ups(:))<0),mean(ups(:)),int8(max(ups(:))<0),max(ups(:)));
    fprintf(1,'         min eps  = %s%4.4f;   mean eps  = %s%4.4f;   max eps  = %s%4.4f;\n'  ,int8(min(eps(:))<0),min(eps(:)),int8(mean(eps(:))<0),mean(eps(:)),int8(max(eps(:))<0),max(eps(:)));
    fprintf(1,'         min tau  = %s%4.4f;   mean tau  = %s%4.4f;   max tau  = %s%4.4f;\n\n',int8(min(tau(:))<0),min(tau(:)),int8(mean(tau(:))<0),mean(tau(:)),int8(max(tau(:))<0),max(tau(:)));

    % update parameters and plot results
    output;
    
    % increment time step
    time = time+dt;
    step = step+1;
    
end

% plot results
if ~isnan(resnorm); output; end

diary off
