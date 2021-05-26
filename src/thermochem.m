% thermochemical solver

if step > 0
    
    % *****  THERMO-CHEMISTRY  ********************************************
    
    % update time step
    Vel    = [U(:)+u(:);W(:)+w(:);fUf(:);fWf(:);fUs(:);fWs(:)];            % combine all velocity components
%     Div_fV = diff(fUs(ic,:),1,2)./h + diff(fWs(:,ic),1,1)./h;            % phase advection/compaction
    Div_fV = diff(fUs(ic,:),1,2)./h + diff(fWs(:,ic),1,1)./h...
            -diff(fUf(ic,:),1,2)./h - diff(fWf(:,ic),1,1)./h;              % phase advection/compaction
    Div_fV = Div_fV/2;
        
    PeTeff = PeT .* (1 - 0.9./(1+exp(-(f(ic,ic)-0.4).*16)));
    PeCeff = PeC .* (1 - 0.9./(1+exp(-(f(ic,ic)-0.4).*16)));

    dt     = CFL*min([ (h/2)^2*max(PeTeff(:)) , h/2/max(abs(Vel)), 0.005./max(abs(Div_fV(:)))]);   % physical time step
    
    % update temperature
    advn_T = flxdiv(T,fUs,fWs,h,ADVN,'adv') + flxdiv(T,fUf,fWf,h,ADVN,'adv'); % advection
        
    cmpt_T = -1/St.*Div_fV;                                                % advection/compaction
        
    diff_T = (diff(T(:,ic),2,1)./h^2 + diff(T(ic,:),2,2)./h^2);            % diffusion
    
    RTin   = (Tin-T)./(5*dt) .* exp(-(L-Z)./(2*h));                        % base injection rate

    RT(ic,ic) = - advn_T - cmpt_T + diff_T./PeTeff + RTin(ic,ic);           % total rate of change
    
    res_T = ((T + f/St)-(To + fo/St))./dt - (theta.*RT + (1-theta).*RTo);  % residual temperature evolution equation
    
    T     = T - res_T.*dt/2;                                               % update temperature solution
    
    if isotherm_topbot                                                     % apply boundary conditions
        T([1 end],:) = To([1 end],:);
    else
        T([1 end],:) = T(ibz,:);
    end
    if isotherm_sides
        T(:,[1 end]) = To(:,[1 end]);
    else
        T(:,[1 end]) = T(:,ibx);
    end
    res_T([1 end],:) = 0;                                                  
    res_T(:,[1 end]) = 0;
    
    
    % update major element composition
    advn_MAJ = flxdiv(MAJs,fUs,fWs,h,ADVN,'adv') + flxdiv(MAJf,fUf,fWf,h,ADVN,'adv'); % advection
        
    cmpt_MAJ = (MAJs(ic,ic)-MAJf(ic,ic)).*Div_fV;                          % advection/compaction
    
    diff_MAJ = (diff(MAJ(:,ic),2,1)./h^2 + diff(MAJ(ic,:),2,2)./h^2);      % diffusion

    RMAJin   = (MAJin-MAJ)./(5*dt) .* exp(-(L-Z)./(2*h));                  % base injection rate

    RMAJ(ic,ic) = - advn_MAJ - cmpt_MAJ + diff_MAJ./PeCeff + RMAJin(ic,ic);    % total rate of change
    
    res_MAJ = (MAJ-MAJo)./dt - (theta.*RMAJ + (1-theta).*RMAJo);           % residual composition evolution equation
    
    MAJ     = MAJ - res_MAJ.*dt/2;                                         % update composition solution
    MAJ     = max(1e-16,min(MAJf,MAJ));                                    % enforce min/max bounds
        
    MAJ([1 end],:) = MAJ(ibz,:);                                           % apply boundary conditions
    MAJ(:,[1 end]) = MAJ(:,ibx);
    res_MAJ([1 end],:) = 0;
    res_MAJ(:,[1 end]) = 0;
    
    
    %*****  PHASE EQUILIBRIUM  ************************************************
    
    % update equilibrium
    [fq,MAJsq,MAJfq]  =  equilibrium(T,MAJ,Pt,perT,perCs,perCf,clap,PhDg);
    
    af = 0.25;
    
    if diseq
        
        % update reaction rate
        RctR_fi = min(Da,0.5/dt).*(fq-f);
        for k  = 1:ceil(kappa)                                             % regularisation
            kk = kappa/ceil(kappa);
            RctR_fi(ic,ic) = RctR_fi(ic,ic) + kk.*(diff(RctR_fi(ic,:),2,2)+diff(RctR_fi(:,ic),2,1))./8;
            RctR_fi([1 end],:) = RctR_fi(ibz,:);
            RctR_fi(:,[1 end]) = RctR_fi(:,ibx);
        end
        RctR_f = af*RctR_fi + (1-af)*RctR_f;

        % update disequilibrium melt fraction
        fin   = (fin-f)./(5*dt) .* exp(-(L-Z)./(2*h));                     % base injection rate
        
        Rf(ic,ic) = + Div_fV + RctR_f(ic,ic) + RMAJin(ic,ic);              % total rate of change
        
        res_f = (f-fo)./dt - (theta.*Rf + (1-theta).*Rfo);                 % residual composition evolution equation
        
        f     = f - res_f.*dt/2;                                           % update composition solution
        f     = max(1e-16,min(1-1e-16,f));                                 % enforce min/max bounds
        
        qind = f<=flim | fq<=flim | (1-f)<=flim | (1-fq)<=flim;
        f(qind)      = fq(qind);                                           % apply equilibrium outside cutoff
        
        f([1 end],:) = f(ibz,:);                                           % apply boundary conditions
        f(:,[1 end]) = f(:,ibx);
        res_f([1 end],:) = 0;
        res_f(:,[1 end]) = 0;
    
    else
        
        % update equilibrium melt fraction
        f              =  af*fq + (1-af)*f;
        
        % update reaction rate
        RctR_fi        =  0.*RctR_f;
        RctR_fi(ic,ic) = (f(ic,ic)-fo(ic,ic))./dt - Div_fV;
        for k  = 1:ceil(kappa)                                             % regularisation
            kk = kappa/ceil(kappa);
            RctR_fi(ic,ic) = RctR_fi(ic,ic) + kk.*(diff(RctR_fi(ic,:),2,2)+diff(RctR_fi(:,ic),2,1))./8;
            RctR_fi([1 end],:) = RctR_fi(ibz,:);
            RctR_fi(:,[1 end]) = RctR_fi(:,ibx);
        end
        RctR_f = af*RctR_fi + (1-af)*RctR_f;
        
    end

    % update phase major element composition
    KMAJ = MAJsq./MAJfq;
    MAJf = MAJ./(f + (1-f).*KMAJ);
    MAJs = MAJ./(f./KMAJ + (1-f));
    MAJs(f<=  flim) = MAJ(f<=  flim);  MAJs(f>=1-flim) = MAJsq(f>=1-flim);
    MAJf(f>=1-flim) = MAJ(f>=1-flim);  MAJf(f<=  flim) = MAJfq(f<=  flim);
        
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
    ISR = double(RctR_f>0).*ISS + double(RctR_f<=0).*ISF;

    
    % *****  TRACE ELEMENTS  **********************************************

    % update incompatible trace element composition
    advn_TRI = flxdiv(TRIs,fUs,fWs,h,ADVN,'adv') + flxdiv(TRIf,fUf,fWf,h,ADVN,'adv'); % advection
        
    cmpt_TRI = (TRIs(ic,ic)-TRIf(ic,ic)).*Div_fV;                          % advection/compaction
    
    diff_TRI = (diff(TRI(:,ic),2,1)./h^2 + diff(TRI(ic,:),2,2)./h^2);      % diffusion
    
    RTRIin   = (TRIin-TRI)./(5*dt) .* exp(-(L-Z)./(2*h));                  % base injection rate

    RTRI(ic,ic) = - advn_TRI - cmpt_TRI + diff_TRI./PeCeff + RTRIin(ic,ic);    % total rate of change
    
    res_TRI = (TRI-TRIo)./dt - (theta.*RTRI + (1-theta).*RTRIo);           % residual composition evolution equation
    
    TRI     = TRI - res_TRI.*dt/2;                                         % update composition solution
    TRI     = max(1e-16,min(TRIf,TRI));                                    % enforce min/max bounds

    TRI([1 end],:) = TRI(ibz,:);                                           % apply boundary conditions
    TRI(:,[1 end]) = TRI(:,ibx);
    res_TRI([1 end],:) = 0;
    res_TRI(:,[1 end]) = 0;
    
    
    % update compatible trace element composition
    advn_TRC = flxdiv(TRCs,fUs,fWs,h,ADVN,'adv') + flxdiv(TRCf,fUf,fWf,h,ADVN,'adv'); % advection
        
    cmpt_TRC = (TRCs(ic,ic)-TRCf(ic,ic)).*Div_fV;                          % advection/compaction
    
    diff_TRC = (diff(TRC(:,ic),2,1)./h^2 + diff(TRC(ic,:),2,2)./h^2);      % diffusion
    
    RTRCin   = (TRCin-TRC)./(5*dt) .* exp(-(L-Z)./(2*h));                  % base injection rate

    RTRC(ic,ic) = - advn_TRC - cmpt_TRC + diff_TRC./PeCeff + RTRCin(ic,ic);    % total rate of change
    
    res_TRC = (TRC-TRCo)./dt - (theta.*RTRC + (1-theta).*RTRCo);           % residual composition evolution equation
    
    TRC     = TRC - res_TRC.*dt/2;                                         % update composition solution
    TRC     = max(1e-16,min(TRCs,TRC));                                    % enforce min/max bounds

    TRC([1 end],:) = TRC(ibz,:);                                           % apply boundary conditions
    TRC(:,[1 end]) = TRC(:,ibx);
    res_TRC([1 end],:) = 0;
    res_TRC(:,[1 end]) = 0;
    

    % *****  RADIOGENIC ISOTOPES  *****************************************

    % update incompatible trace element composition
    advn_IRP = flxdiv(IRPs,fUs,fWs,h,ADVN,'adv') + flxdiv(IRPf,fUf,fWf,h,ADVN,'adv'); % advection
        
    cmpt_IRP = (IRPs(ic,ic)-IRPf(ic,ic)).*Div_fV;                          % advection/compaction
    
    diff_IRP = (diff(IRP(:,ic),2,1)./h^2 + diff(IRP(ic,:),2,2)./h^2);      % diffusion
    
    RIRPin   = (IRPin.*exp(-time./DIRP.*log(2))-IRP)./(5*dt) .* exp(-(L-Z)./(2*h)); % base injection rate
        
    RIRP(ic,ic) = - advn_IRP - cmpt_IRP + diff_IRP./PeCeff - DcyR_IR(ic,ic) + RIRPin(ic,ic); % total rate of change;
    
    res_IRP = (IRP-IRPo)./dt - (theta.*RIRP + (1-theta).*RIRPo);           % residual composition evolution equation
    
    IRP     = IRP - res_IRP.*dt/2;                                         % update composition solution
    IRP     = max(1e-16,min(1e3,IRP));                                     % enforce min/max bounds

    IRP([1 end],:) = IRP(ibz,:);                                           % apply boundary conditions
    IRP(:,[1 end]) = IRP(:,ibx);
    res_IRP([1 end],:) = 0;
    res_IRP(:,[1 end]) = 0;
    
    
    % update compatible trace element composition
    advn_IRD = flxdiv(IRDs,fUs,fWs,h,ADVN,'adv') + flxdiv(IRDf,fUf,fWf,h,ADVN,'adv'); % advection
        
    cmpt_IRD = (IRDs(ic,ic)-IRDf(ic,ic)).*Div_fV;                          % advection/compaction
    
    diff_IRD = (diff(IRD(:,ic),2,1)./h^2 + diff(IRD(ic,:),2,2)./h^2);      % diffusion
    
    RIRDin   = (IRDin+IRPin.*(1-exp(-time./DIRP.*log(2)))-IRD)./(5*dt) .* exp(-(L-Z)./(2*h)); % base injection rate

    RIRD(ic,ic) = - advn_IRD - cmpt_IRD + diff_IRD./PeCeff + DcyR_IR(ic,ic) + RIRDin(ic,ic); % total rate of change;
    
    res_IRD = (IRD-IRDo)./dt - (theta.*RIRD + (1-theta).*RIRDo);           % residual composition evolution equation
    
    IRD     = IRD - res_IRD.*dt/2;                                         % update composition solution
    IRD     = max(1e-16,min(1e3,IRD));                                     % enforce min/max bounds

    IRD([1 end],:) = IRD(ibz,:);                                           % apply boundary conditions
    IRD(:,[1 end]) = IRD(:,ibx);
    res_IRD([1 end],:) = 0;
    res_IRD(:,[1 end]) = 0;
    
    
    % *****  STABLE ISOTOPES  *********************************************

    % update solid stable isotope composition
    advn_ISS = flxdiv(ISS,U,W,h,ADVN,'adv');                               % advection/compaction
    
    diff_ISS = (diff(ISS(:,ic),2,1)./h^2 + diff(ISS(ic,:),2,2)./h^2);      % diffusion
    
    RctR_ISS = -(ISR-ISS).*RctR_f./max(flim,1-f);                          % reactive transfer rate

    RISSin   = (ISSin-ISS)./(5*dt) .* exp(-(L-Z)./(2*h));                  % base injection rate

    RISS(ic,ic) = - advn_ISS + diff_ISS./PeCeff + RctR_ISS(ic,ic) + RISSin(ic,ic);  % total rate of change
    
    res_ISS = (ISS-ISSo)./dt - (theta.*RISS + (1-theta).*RISSo);           % residual composition evolution equation
    
    ISS     = ISS - res_ISS.*dt/2;                                         % update composition solution
    ISS     = max(min(ISSin(:)),min(max(ISSin(:)),ISS));                   % enforce min/max bounds

    ISS([1 end],:) = ISS(ibz,:);                                           % apply boundary conditions
    ISS(:,[1 end]) = ISS(:,ibx);    
    res_ISS([1 end],:) = 0;
    res_ISS(:,[1 end]) = 0;
    
    
    % update fluid stable isotope composition
    advn_ISF = flxdiv(ISF,Uf,Wf,h,ADVN,'adv');                             % advection/compaction
    
    diff_ISF = (diff(ISF(:,ic),2,1)./h^2 + diff(ISF(ic,:),2,2)./h^2);      % diffusion
    
    RctR_ISF = (ISR-ISF).*RctR_f./max(flim,f);                             % reactive transfer rate
    
    RISFin   = (ISFin-ISF)./(5*dt) .* exp(-(L-Z)./(2*h));                  % base injection rate

    RISF(ic,ic) = - advn_ISF + diff_ISF./PeCeff + RctR_ISF(ic,ic) + RISFin(ic,ic); % total rate of change
    
    res_ISF = (ISF-ISFo)./dt - (theta.*RISF + (1-theta).*RISFo);           % residual composition evolution equation
    
    ISF     = ISF - res_ISF.*dt/2;                                         % update composition solution
    ISF     = max(min(ISFin(:)),min(max(ISFin(:)),ISF));                   % enforce min/max bounds

    ISF([1 end],:) = ISF(ibz,:);                                           % apply boundary conditions
    ISF(:,[1 end]) = ISF(:,ibx);    
    res_ISF([1 end],:) = 0;
    res_ISF(:,[1 end]) = 0;
    
    
    %*****  DAMAGE EVOLUTION  *************************************************
    
    epsDMG   = max(1e-16,eps-tau./10.^etav);                               % failure damage rate

    advn_DMG = flxdiv(DMG,U,W,h,ADVN,'adv');                               % damage advection
    
    diff_DMG = (diff(DMG(:,ic),2,1)./h^2 + diff(DMG(ic,:),2,2)./h^2);      % regularisation
        
    RDMG(ic,ic) = -advn_DMG + diff_DMG./PeCeff ...
                + epsDMG(ic,ic) - RHEAL.*DMG(ic,ic);                       % total rate of change
    
    res_DMG = (DMG-DMGo)./dt - (theta.*RDMG + RDMGo);                      % residual damage evolution equation
        
    DMG = DMG - res_DMG.*dt/2;                                               % update composition solution
    DMG = max(1e-16,DMG);                                                  % enforce min bound
    
    DMG([1 end],:) = DMG(ibz,:);                                           % apply boundary conditions
    DMG(:,[1 end]) = DMG(:,ibx);
    res_DMG([1 end],:) = 0;
    res_DMG(:,[1 end]) = 0;
    
else
    
    % update equilibrium
    [fq,MAJsq,MAJfq]  =  equilibrium(T,MAJ,Pt,perT,perCs,perCf,clap,PhDg);
    f                 = fq;
    
end
