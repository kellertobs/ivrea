%*****  TWO-PHASE RHEOLOGY  ***********************************************

% update tensor magnitudes
eps(ic,ic) = (  (exx(ic,ic).^2 + ezz(ic,ic).^2 ...
           + 2.*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2).*0.25)./2).^0.5 + 1e-16;
eps([1 end],:) = eps(ibz,:);                                               % periodic boundaries
eps(:,[1 end]) = eps(:,ibx);

tau(ic,ic) = (  (txx(ic,ic).^2 + tzz(ic,ic).^2 ...
           + 2.*(txz(1:end-1,1:end-1).^2 + txz(2:end,1:end-1).^2 + txz(1:end-1,2:end).^2 + txz(2:end,2:end).^2).*0.25)./2).^0.5 + 1e-16;
tau([1 end],:) = tau(ibz,:);                                               % periodic boundaries
tau(:,[1 end]) = tau(:,ibx);


plim = double(f>=flim);
for k  = 1:ceil(kappa)                                                     % regularisation
    kk = kappa/ceil(kappa);
    plim(ic,ic) = plim(ic,ic) + kk.*(diff(plim(ic,:),2,2)+diff(plim(:,ic),2,1))./8;
    plim([1 end],:) = plim(ibz,:);
    plim(:,[1 end]) = plim(:,ibx);
end

Pe = p.*plim + Pt.*(1-plim);
yieldt_GM = max(1e-16, 1 + Pe  );
yieldt_MC = max(1e-16, 2 + Pe/2);
yieldt    = min(yieldt_GM,yieldt_MC) .* 0.25.^DMG + etamin.*eps + bnchmrk*10;

% update rheological parameters
etav  = exp(Es*(1./T-1./T0) -30.*(f-f0)) .* (1/2+1/2*(eps./eps0).^-n); %.* 0.01.^(MAJ-MAJ0)
etav  = (1./etamax + 1./etav).^-1 + etamin;
etav  = log10(etav);
etav([1 end],:) = etav(ibz,:);
etav(:,[1 end]) = etav(:,ibx);
    
etay  =  log10(yieldt)-log10(eps);                                         % shear visco-plasticity
etay  =  min(etav,etay);

limf  =  max(flim,min(1-flim,f));
zetay =  etay - log10(limf.*(1-limf).^0.5);

for k  = 1:ceil(kappa)                                                     % regularisation
    kk = kappa/ceil(kappa);
    etay(ic,ic) = etay(ic,ic) + kk.*(diff(etay(ic,:),2,2)+diff(etay(:,ic),2,1))./8;
    etay([1 end],:) = etay(ibz,:);
    etay(:,[1 end]) = etay(:,ibx);
end

for k  = 1:ceil(kappa)                                                     % regularisation
    kk = kappa/ceil(kappa);
    zetay(ic,ic) = zetay(ic,ic) + kk.*(diff(zetay(ic,:),2,2)+diff(zetay(:,ic),2,1))./8;
    zetay([1 end],:) = zetay(ibz,:);
    zetay(:,[1 end]) = zetay(:,ibx);
end

 eta  =   etay.*(1-gamma) +  eta.*gamma;                                            % effective shear viscosity
zeta  =  zetay.*(1-gamma) + zeta.*gamma;                                            % effective shear viscosity

etac  = (eta(im,im)+eta(ip,im) + eta(im,ip)+eta(ip,ip)).*0.25;             % interpolate to cell corners

K   = log10((limf/f0).^3 .* (1-limf).^1.5 .* exp(-Ef*(1./T-1./T0))) + KDMG.*DMG;  % segregation coefficient

% update iterative and physical time step sizes
dtW = (10.^(( eta(im,:)+ eta(ip,:)).*0.5)./(h/2)^2 ...
    +  10.^((zeta(im,:)+zeta(ip,:)).*0.5)./(h/2)^2).^-1;                   % W iterative step size
dtU = (10.^(( eta(:,im)+ eta(:,ip)).*0.5)./(h/2)^2 ...
    +  10.^((zeta(:,im)+zeta(:,ip)).*0.5)./(h/2)^2).^-1;                   % U iterative step size
dtP = 1/4*(1./(10.^eta) + 10.^K./(h/2)^2).^-1;                             % P iterative step size


if step > 0
    
    %*****  DAMAGE EVOLUTION  *************************************************
    
    epsDMG   = max(1e-16,eps-tau./10.^etav);                               % failure damage rate

    V_GrdDMG = flxdiv(DMG,U+UBG,W+WBG,h,'FROMM','adv');                    % damage advection
    
    Lpl_DMG  = (diff(DMG(:,ic),2,1)./h^2 + diff(DMG(ic,:),2,2)./h^2);      % regularisation
        
    RDMG(ic,ic) = -V_GrdDMG + Lpl_DMG/PeC ...
                + epsDMG(ic,ic) - RHEAL.*DMG(ic,ic);                       % total rate of change
    
    res_DMG = (DMG-DMGo)./dt - (theta.*RDMG + RDMGo);                      % residual damage evolution equation
        
    DMG = DMG - res_DMG.*dt/10;                                            % update composition solution
    
    DMG([1 end],:) = DMG(ibz,:);                                           % apply boundary conditions
    DMG(:,[1 end]) = DMG(:,ibx);
    res_DMG([1 end],:) = 0;
    res_DMG(:,[1 end]) = 0;
    
    
    %*****  PHASE EQUILIBRIUM  ************************************************
        
    % update equilibrium
    [fq,MAJsq,MAJfq]  =  equilibrium(T,MAJ,perT,perCs,perCf,PhDg);
    RctR_fi = min(Da,1/(2*dt)).*(fq-f);
    for k  = 1:ceil(kappa)                                                 % regularisation
        kk = kappa/ceil(kappa);
        RctR_fi(ic,ic) = RctR_fi(ic,ic) + kk.*(diff(RctR_fi(ic,:),2,2)+diff(RctR_fi(:,ic),2,1))./8;
        RctR_fi([1 end],:) = RctR_fi(ibz,:);
        RctR_fi(:,[1 end]) = RctR_fi(:,ibx);
    end
    RctR_f = (1-gamma).*RctR_fi + gamma*RctR_f;

    % update phase major element composition
    %KMAJq = MAJsq./MAJfq;
    MAJf = MAJfq;%MAJ./(f + (1-f).*KMAJq);
    MAJs = MAJsq;%MAJ./(f./KMAJq + (1-f));
    MAJf(f<=  flim) = MAJfq(f<=  flim);  MAJs(f<=  flim) = MAJ(f<=    flim);
    MAJf(f>=1-flim) = MAJ  (f>=1-flim);  MAJs(f>=1-flim) = MAJsq(f>=1-flim);
    
    % update incompatible trace element phase compositions
    TRIf = TRI./(f + (1-f).*KTRI);
    TRIs = TRI./(f./KTRI + (1-f));
    TRIs(f<=  flim) = TRI(f<=  flim);
    TRIf(f>=1-flim) = TRI(f>=1-flim);
    
    % update compatible trace element phase compositions
    TRCf = TRC./(f + (1-f).*KTRC);
    TRCs = TRC./(f./KTRC + (1-f));
    TRCs(f<=  flim) = TRC(f<=  flim);
    TRCf(f>=1-flim) = TRC(f>=1-flim);

    % update radiogenic parent isotope phase compositions
    IRPf = IRP./(f + (1-f).*KIRP);
    IRPs = IRP./(f./KIRP + (1-f));
    IRPs(f<=  flim) = IRP(f<=  flim);
    IRPf(f>=1-flim) = IRP(f>=1-flim);
    
    % update radiogenic daughter isotope phase compositions
    IRDf = IRD./(f + (1-f).*KIRD);
    IRDs = IRD./(f./KIRD + (1-f));
    IRDs(f<=  flim) = IRD(f<=  flim);
    IRDf(f>=1-flim) = IRD(f>=1-flim);
    
    % radiogenic isotope decay rate
    DcyR_IR = IRP./DIRP.*log(2);
    
    % stable isotope transfer composition
    ISR = double(RctR_f>=0).*ISS + double(RctR_f<0).*ISF;
    
end