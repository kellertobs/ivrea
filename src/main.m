% print run header
fprintf('\n\n*****  ivrea  |  %s  |  %s  *****\n\n',runID,datetime);

% load custom colormap
load ocean;  

ic  =  2:N-1 ;
ip  =  2:N   ;
im  =  1:N-1 ;
ibz = [2 N-1];  % boundary indices for closed boundary condition
ibx = [N-1 2];  % boundary indices for periodic boundary condition

% produce smooth random perturbations
rng(5);
dr = randn(N,N);
for i = 1:max(smx,smz)
    dr(ic,ic) = dr(ic,ic) ...
              + smz./max(smx,smz).*diff(dr(:,ic),2,1)./8 ...
              + smx./max(smx,smz).*diff(dr(ic,:),2,2)./8;
    dr = dr - mean(dr(:)); 
    dr = dr./max(abs(dr(:)));
    dr([1 end],:) = dr([end-1 2],:);
    dr(:,[1 end]) = dr(:,[end-1 2]);
end

% get coordinate arrays
z     = -L/2-h/2:h:L/2+h/2;
x     = -L/2-h/2:h:L/2+h/2;
xc    = (x(im)+x(ip))./2;
zc    = (z(im)+z(ip))./2;
[X,Z] = meshgrid(x,z);
Xbnd = X; Zbnd = Z;

% set gaussian perturbation shape function
gs   = exp(-(X+xpos).^2./wx^2).*exp(-(Z+zpos).^2./wz^2);

% initialise background velocity fields
Pu  = Pu+1e-16;  Si = Si+1e-16;
WP  =  (Z(im,:)+Z(ip,:))/2/L*Pu*L;  % pure   shear WBG
WS  = -(X(im,:)+X(ip,:))/2/L*Si*L;  % simple shear WBG
UP  = -(X(:,im)+X(:,ip))/2/L*Pu*L;  % pure   shear UBG
US  = -(Z(:,im)+Z(:,ip))/2/L*Si*L;  % simple shear UBG
WS  = 0.*WS; US = 2.*US; 
WBG = WP+WS;
UBG = UP+US;     
clear  WP WS UP US

% initialise thermo-chemical solution and residual fields
if bnchmrk
    mms;
    f   = f_mms;
    T   = T_mms;
    MAJ = C_mms;
else
    T         =  Tc-(T0-Tc).*erf((-Z-L/2)./D).*(1 + T1.*dr + T2.*gs); Tin = T;
    MAJ       =  MAJ0 .* (1 + MAJ1.*dr + MAJ2.*gs); MAJin = MAJ;
    TRC       =  TRC0 .* (1 + TRC1.*dr + TRC2.*gs); TRCin = TRC;
    SIS       =  SIS0 .* (1 + SIS1.*dr + SIS2.*gs); SISin = SIS;
    [f,MAJs,MAJf] =  equilibrium(T ,MAJ ,perT,perCs,perCf,PhDg);  fin = f;
    f0        =  equilibrium(T0,MAJ0,perT,perCs,perCf,PhDg);
end
fprintf('\n\n*****  initial condition: T0 = %1.3f;  C0 = %1.3f;  f0 = %1.3f;\n\n',T0,MAJ0,f0);
res_f   = 0.*f; df = 0.*f; Lpl_f = 0.*f;
res_T   = 0.*T; dT = 0.*T;
res_MAJ = 0.*MAJ; dMAJ = 0.*MAJ;
res_TRC = 0.*TRC;
res_SIS = 0.*SIS;

% initialise mechanical solution and residual fields
W     =  0.*WBG;  res_W = 0.*W;  dW = 0.*W;
U     =  0.*UBG;  res_U = 0.*U;  dU = 0.*U;
P     =  0.*f;    res_P = 0.*P;  dP = 0.*P;
u     =  0.*U;    uf = U+u./((f(:,im)+f(:,ip))./2);
w     =  0.*W;    wf = W+w./((f(im,:)+f(ip,:))./2);
p     =  0.*P;

% initialise parameter fields
V_GrdT    =  0.*T;    Lpl_T = 0.*T;
Div_fMAJV =  0.*MAJ;  Lpl_MAJ = 0.*MAJ;
Div_fTRCV =  0.*TRC;  Lpl_TRC = 0.*TRC;
Div_fSISV =  0.*SIS;  Lpl_SIS = 0.*SIS;
RctR_f    =  0.*f;    RctR_fo = RctR_f;
ups       =  0.*P;  upss = 0.*P;  Div_fV = 0.*P;  Div_fVBG = 0.*P;
eps0      =  max(B/L,abs(Pu) + abs(Si)) + 1e-16;  
exx       =  0.*P - Pu;  ezz = 0.*P + Pu;  exz = zeros(N-1,N-1) - Si;  eps = 0.*P + (abs(Pu) + abs(Si));  
txx       =  0.*exx;  tzz = 0.*ezz;  txz = 0.*exz;  tau = 0.*eps;
eta       =  exp(Es*(1./T-1./T0)) .* (1-(f-f0)./0.3).^8;
eta       =  (1./etamax + 1./eta ).^-1 + etamin;
eta       =  log10( eta );
zeta      =  eta - max(-6,log10(f));
yieldt    =  ones(size(P));
rctr      =  ones(size(MAJ));

% initialise timing parameters
step  =  0;
time  =  0;
dt    =  1e-16;

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
    To   = T;   V_GrdTo    = V_GrdT;    Lpl_To   = Lpl_T;
    MAJo = MAJ; Div_fMAJVo = Div_fMAJV; Lpl_MAJo = Lpl_MAJ;
    TRCo = TRC; Div_fTRCVo = Div_fTRCV; Lpl_TRCo = Lpl_TRC;
    SISo = SIS; Div_fSISVo = Div_fSISV; Lpl_SISo = Lpl_SIS;
    fo   = f;   Div_fVo    = Div_fV;    Lpl_fo   = Lpl_f;    RctR_fo = RctR_f;
    dto   = dt;
    
    % reset residual norms and iteration count
    resnorm  = 1e3;
    resnorm0 = resnorm;
    it       = 0;
    
    % initialise iterative solution guess
    dWi = 0.*W;
    dUi = 0.*U;
    dPi = 0.*P;
    
    % non-linear iteration loop
    while resnorm/resnorm0 >= rtol && resnorm >= atol && it <= maxit || it <= minit
                
        % store previous iterative solution guess  
        Wi = W;
        Ui = U;
        Pi = P;
        
        % update thermo-chemical evolution
        if ~mod(it,nup)
            
            up2date;
            
            Div_fV(ic,ic) = flxdiv(1-f,U+UBG,W+WBG,h,'flx');               % phase advection/compaction

            if step > 0
                
            % update melt fraction
            Lpl_f(ic,ic)  = (diff(f(:,ic),2,1)./h^2 ...
                          +  diff(f(ic,:),2,2)./h^2);                      % diffusion
                      
            RbndF  = 0;%max(0,-(f-fin)/5/dt .* exp(-(( Z)-L/2).^2./(L/20)^2)) ...
                   %+ min(0, (f-fin)/5/dt .* exp(-((-Z)-L/2).^2./(L/20)^2));
            
            res_f = (f-fo)./dt - (  theta .*(Div_fV  + Lpl_f /PeC + RctR_f ) ...
                               + (1-theta).*(Div_fVo + Lpl_fo/PeC + RctR_fo) ...
                               + RbndF);                                   % residual phase evolution equation
                           
            df = -res_f.*dt/8;
            
            df([1 end],:)  = 0;                                        % apply boundary conditions
            df(:,[1 end])  = df(:,ibx);
                    
            f = f + alpha*df;                                           % update phase solution
            f = max(1e-16,min(1-1e-16,f));
            f([1 end],:) = fq([1 end],:);
            
            % update temperature
            V_GrdT(ic,ic) = flxdiv(T,U+u+UBG,W+w+WBG,h,'adv');             % advection
            
            Lpl_T(ic,ic)  = (diff(T(:,ic),2,1)./h^2 ...
                          +  diff(T(ic,:),2,2)./h^2);                      % diffusion
                                
            RbndT  = 0;%max(0,-(T-Tin)/5/dt .* exp(-(( Z)-L/2).^2./(L/20)^2)) ...
                   %+ min(0,-(T-Tin)/5/dt .* exp(-((-Z)-L/2).^2./(L/20)^2));
                
            res_T = (T-To)./dt - (  theta .*(-V_GrdT  + Lpl_T /PeT - RctR_f /St) ...
                               + (1-theta).*(-V_GrdTo + Lpl_To/PeT - RctR_fo/St) ...
                               + RbndT);                                   % residual temperature evolution equation
            
            dT = -res_T.*dt/8;
            
            dT([1 end],:)  = dT(ibz,:);                                % apply boundary conditions
            dT(:,[1 end])  = dT(:,ibx);
                             
            T = T + alpha*dT;                                           % update temperature solution
            
            % update major element composition
            Div_fMAJV(ic,ic) = flxdiv((1-f).*MAJs,U +UBG,W +WBG,h,'flx') ...
                             + flxdiv(   f .*MAJf,uf+UBG,wf+WBG,h,'flx');  % advection/compaction
            
            Lpl_MAJ(ic,ic)  = (diff(MAJ(:,ic),2,1)./h^2 ...
                            +  diff(MAJ(ic,:),2,2)./h^2);                  % diffusion

            RbndMAJ = 0;%max(0,-(MAJ-MAJin)/5/dt .* exp(-(( Z)-L/2).^2./(L/20)^2)) ...
                    %+ min(0,-(MAJ-MAJin)/5/dt .* exp(-((-Z)-L/2).^2./(L/20)^2));

            res_MAJ = (MAJ-MAJo)./dt - (  theta .*(-Div_fMAJV  + Lpl_MAJ /PeC) ...
                                     + (1-theta).*(-Div_fMAJVo + Lpl_MAJo/PeC) ...
                                     + RbndMAJ);                           % residual composition evolution equation
            
            dMAJ = -res_MAJ.*dt/8;
            
            dMAJ([1 end],:)  = dMAJ(ibz,:);                                % apply boundary conditions
            dMAJ(:,[1 end])  = dMAJ(:,ibx);
                        
            MAJ = MAJ + alpha*dMAJ;                                     % update composition solution
            MAJ = max(1e-16,min(1-1e-16,MAJ));

            % update trace element composition
            Div_fTRCV(ic,ic) = flxdiv((1-f).*TRCs,U +UBG,W +WBG,h,'flx') ...
                             + flxdiv(   f .*TRCf,uf+UBG,wf+WBG,h,'flx');  % advection/compaction
            
            Lpl_TRC(ic,ic)  = (diff(TRC(:,ic),2,1)./h^2 ...
                            +  diff(TRC(ic,:),2,2)./h^2);                  % diffusion
                        
            RbndTRC = 0;%max(0,-(TRC-TRCin)/5/dt .* exp(-(( Z)-L/2).^2./(L/20)^2)) ...
                    %+ min(0,-(TRC-TRCin)/5/dt .* exp(-((-Z)-L/2).^2./(L/20)^2));
            
            res_TRC = (TRC-TRCo)./dt - (  theta .*(-Div_fTRCV  + Lpl_TRC /PeC) ...
                                     + (1-theta).*(-Div_fTRCVo + Lpl_TRCo/PeC) ...
                                     + RbndTRC);                           % residual composition evolution equation
            
            res_TRC([1 end],:) = res_TRC(ibz,:);                           % apply boundary conditions
            res_TRC(:,[1 end]) = res_TRC(:,ibx);
                        
            TRC = TRC - res_TRC.*dt/8;                                     % update composition solution
            
            % update stable isotope composition
            Div_fSISV(ic,ic) = flxdiv((1-f).*SIS,U +UBG,W +WBG,h,'flx') ...
                             + flxdiv(   f .*SIS,uf+UBG,wf+WBG,h,'flx');   % advection/compaction
            
            Lpl_SIS(ic,ic)  = (diff(SIS(:,ic),2,1)./h^2 ...
                            +  diff(SIS(ic,:),2,2)./h^2);                  % diffusion
                        
            RbndSIS = 0;%max(0,-(SIS-SISin)/5/dt .* exp(-(( Z)-L/2).^2./(L/20)^2)) ...
                    %+ min(0,-(SIS-SISin)/5/dt .* exp(-((-Z)-L/2).^2./(L/20)^2));
            
            res_SIS = (SIS-SISo)./dt - (  theta .*(-Div_fSISV  + Lpl_SIS /PeC) ...
                                     + (1-theta).*(-Div_fSISVo + Lpl_SISo/PeC) ...
                                     + RbndSIS);                           % residual composition evolution equation
            
            res_SIS([1 end],:) = res_SIS(ibz,:);                           % apply boundary conditions
            res_SIS(:,[1 end]) = res_SIS(:,ibx);
                        
            SIS = SIS - res_SIS.*dt/8;                                     % update composition solution
            end
        end
        
        % update segregation velocities and compaction pressure
        w   = -(K(im,:)+K(ip,:)).*0.5 .* (diff(P,1,1)./h + B);             % segregation z-velocity
        
        w(:,[1 end]) = w(:,ibx);
                
        u   = -(K(:,im)+K(:,ip)).*0.5 .* (diff(P,1,2)./h);                 % segregation x-velocity
        
        u([1 end],:) = u(ibz,:);                                       % apply boundary conditions
        u(:,[1 end]) = [sum(u(:,[1 end]),2)./2, ...
                        sum(u(:,[1 end]),2)./2];
                
        uf = U+u./((f(:,im)+f(:,ip))./2);                                  % get fluid x-velocity      
        wf = W+w./((f(im,:)+f(ip,:))./2);                                  % get fluid z-velocity
            
        p  = -10.^zeta .* ups;                                             % compaction pressure
        
        p([1 end],:) = 0;                                                  % apply boundary conditions
        p(:,[1 end]) = p(:,ibx);
                
        % update strain rates
        exx(:,ic)   = diff(U,1,2)./h - ups(:,ic)./3 - Pu;                  % x-normal strain rate
        exx([1 end],:)   = exx(ibz,:);                                     % apply boundary conditions
        exx(:,[1 end])   = exx(:,ibx);               
        ezz(ic,:)   = diff(W,1,1)./h - ups(ic,:)./3 + Pu;                  % z-normal strain rate
        ezz([1 end],:)   = ezz(ibz,:);                                     % apply boundary conditions
        ezz(:,[1 end])   = ezz(:,ibx);          
        exz              = 1/2.*(diff(U,1,1)./h+diff(W,1,2)./h) - Si;      % shear strain rate
        
        % update stresses
        txx = 10.^eta .* exx;                                              % x-normal stress
        tzz = 10.^eta .* ezz;                                              % z-normal stress
        txz = 10.^etac.* exz;                                              % xz-shear stress  
        
        % update z-reference velocity
        Div_tz = diff(tzz(:,ic),1,1)./h + diff(txz,1,2)./h;                % z-stress divergence
        
        res_W(:,ic) = - Div_tz + diff(P(:,ic),1,1)./h ...                  % residual z-momentum equation
                               + diff(p(:,ic),1,1)./h; % - ((f(im,ic)+f(ip,ic))./2-f0).*B;
                                
        if bnchmrk; res_W = res_W - src_W_mms; end
        
        dW = -res_W.*dtW;
                
        dW(:,[1 end]) = dW(:,ibx);
        
        W = Wi + alpha*dW + beta*dWi;                                      % update z-velocity solution

        dWi = W - Wi;
        
        % update x-reference velocity        
        Div_tx  = diff(txx(ic,:),1,2)./h + diff(txz,1,1)./h;               % x-stress divergence
        
        res_U(ic,:) = - Div_tx + diff(P(ic,:),1,2)./h ...                  % residual x-momentum equation
                                    + diff(p(ic,:),1,2)./h;
        if bnchmrk; res_U = res_U - src_U_mms; end
        
        dU = -res_U.*dtU;
        
        dU([1 end],:) = dU(ibz,:);                                         % apply boundary conditions
        dU(:,[1 end]) = [sum(dU(:,[1 end]),2)./2, ...
                         sum(dU(:,[1 end]),2)./2];

        dU = dU-mean(dU(:));                                               % remove mean from update
        
        U = Ui + alpha*dU + beta*dUi;                                      % update x-velocity solution
      
        dUi = U - Ui;

        % update velocity divergences
        ups(ic,ic) = diff(U(ic,:),1,2)./h + diff(W(:,ic),1,1)./h;          % matrix velocity divergence
                             
        ups([1 end],:) = ups(ibz,:);                                       % apply boundary conditions
        ups(:,[1 end]) = ups(:,ibx);
        
        upss(ic,ic) = diff(u(ic,:),1,2)./h + diff(w(:,ic),1,1)./h;         % segregation velocity divergence
        upss([1 end],:) = upss(ibz,:);                                     % apply boundary conditions    
        upss(:,[1 end]) = upss(:,ibx);          
        
        % update reference pressure
        res_P = ups + upss;                                                % residual mass equation
        if bnchmrk; res_P = res_P - src_P_mms; end
        
        dP = -res_P.*dtP;
        
        dP([1 end],:) = 0;                                                 % apply boundary conditions
        dP(:,[1 end]) = dP(:,ibx);
        
        dP([1 end],:) = 0;                                             % apply boundary conditions
        dP(:,[1 end]) = dP(:,ibx);
        
        dP = dP-mean(dP(:));                                               % remove mean from update
            
        P = Pi + alpha*dP + beta*dPi;                                      % update pressure solution

        dPi = P - Pi;
        
        % check and report convergence every max(100,nup) iterations
        if ~mod(it,max(100,nup)); report; end
        
        it = it+1;

    end
    
    % clean workspace
    clear Wi Ui Pi dtWi dtUi dtPi dtWii dtUii dtPii fo Div_fVo Div_fVBGo Div_tz Div_tx dtW dtU dtP etac
    
    % print diagnostics
    fprintf(1,'\n         time to solution = %4.4f sec\n\n',toc);

    fprintf(1,'         min U    = %s%4.4f;   mean U    = %s%4.4f;   max U    = %s%4.4f;\n'  ,int8(min(  U(:))<0),min(  U(:)),int8(mean(  U(:))<0),mean(  U(:)),int8(max(  U(:))<0),max(  U(:)));
    fprintf(1,'         min W    = %s%4.4f;   mean W    = %s%4.4f;   max W    = %s%4.4f;\n'  ,int8(min( -W(:))<0),min( -W(:)),int8(mean( -W(:))<0),mean( -W(:)),int8(max( -W(:))<0),max( -W(:)));
    fprintf(1,'         min P    = %s%4.4f;   mean P    = %s%4.4f;   max P    = %s%4.4f;\n\n',int8(min(  P(:))<0),min(  P(:)),int8(mean(  P(:))<0),mean(  P(:)),int8(max(  P(:))<0),max(  P(:)));

    fprintf(1,'         min u    = %s%4.4f;   mean u    = %s%4.4f;   max u    = %s%4.4f;\n'  ,int8(min(  u(:))<0),min(  u(:)),int8(mean(  u(:))<0),mean(  u(:)),int8(max(  f(:))<0),max(  u(:)));
    fprintf(1,'         min w    = %s%4.4f;   mean w    = %s%4.4f;   max w    = %s%4.4f;\n'  ,int8(min( -w(:))<0),min( -w(:)),int8(mean( -w(:))<0),mean( -w(:)),int8(max(  f(:))<0),max( -w(:)));
    fprintf(1,'         min p    = %s%4.4f;   mean p    = %s%4.4f;   max p    = %s%4.4f;\n\n',int8(min(  p(:))<0),min(  p(:)),int8(mean(  p(:))<0),mean(  p(:)),int8(max(  f(:))<0),max(  p(:)));
    
    fprintf(1,'         min f    = %s%4.4f;   mean f    = %s%4.4f;   max f    = %s%4.4f;\n'  ,int8(min(  f(:))<0),min(  f(:)),int8(mean(  f(:))<0),mean(  f(:)),int8(max(  f(:))<0),max(  f(:)));
    fprintf(1,'         min K    = %s%4.4f;   mean K    = %s%4.4f;   max K    = %s%4.4f;\n\n',int8(min(  K(:))<0),min(  K(:)),int8(mean(  K(:))<0),mean(  K(:)),int8(max(  K(:))<0),max(  K(:)));
    
    fprintf(1,'         min  eta = %s%4.4f;   mean  eta = %s%4.4f;   max  eta = %s%4.4f;\n'  ,int8(min( eta(:))<0),min( eta(:)),int8(mean( eta(:))<0),mean( eta(:)),int8(max( eta(:))<0),max( eta(:)));
    fprintf(1,'         min zeta = %s%4.4f;   mean zeta = %s%4.4f;   max zeta = %s%4.4f;\n\n',int8(min(zeta(:))<0),min(zeta(:)),int8(mean(zeta(:))<0),mean(zeta(:)),int8(max(zeta(:))<0),max(zeta(:)));
    
    fprintf(1,'         min ups  = %s%4.4f;   mean ups  = %s%4.4f;   max ups  = %s%4.4f;\n'  ,int8(min(ups(:))<0),min(ups(:)),int8(mean(ups(:))<0),mean(ups(:)),int8(max(ups(:))<0),max(ups(:)));
    fprintf(1,'         min eps  = %s%4.4f;   mean eps  = %s%4.4f;   max eps  = %s%4.4f;\n'  ,int8(min(eps(:))<0),min(eps(:)),int8(mean(eps(:))<0),mean(eps(:)),int8(max(eps(:))<0),max(eps(:)));
    fprintf(1,'         min tau  = %s%4.4f;   mean tau  = %s%4.4f;   max tau  = %s%4.4f;\n\n',int8(min(tau(:))<0),min(tau(:)),int8(mean(tau(:))<0),mean(tau(:)),int8(max(tau(:))<0),max(tau(:)));

    % update parameters and plot results
    if ~mod(step,nop); up2date; output; end
    
    % increment time step
    time = time+dt;
    step = step+1;
    
end

% plot results
if ~isnan(resnorm); up2date; output; end

diary off
