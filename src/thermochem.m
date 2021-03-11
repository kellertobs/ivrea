% thermochemical solver

if step > 0

    % *****  THERMO-CHEMISTRY  ********************************************
    
    % update melt fraction
    Div_fV = flxdiv(1-f,U+UBG,W+WBG,h,'FLXDIV','flx');                     % phase advection/compaction
    
    Vel    = [U(:)+UBG(:);W(:)+WBG(:);u(:);w(:);uf(:);wf(:)];              % combine all velocity components
    dt     = CFL*min([h/2/max(abs(Vel)), 0.0025./max(abs(Div_fV(:)))]);    % physical time step

    Lpl_f  = (diff(f(:,ic),2,1)./h^2 + diff(f(ic,:),2,2)./h^2);            % diffusion
    
    Rf(ic,ic) =  Div_fV + Lpl_f/PeC + RctR_f(ic,ic);                       % total rate of change
    
    res_f = (f-fo)./dt - (theta.*Rf + (1-theta).*Rfo);                     % residual phase evolution equation
    
    f     = f - res_f.*dt/10;                                              % update phase fraction
    
    iq    = fq<=flim | fq>=(1-flim);
    f(iq) = fq(iq);  res_f(iq) = 0;                                        % enforce equilibrium at low
    

    f([1 2 end-1 end],:) = fq([1 2 end-1 end],:);                          % apply boundary conditions
    if abs(Si) > abs(Pu)
        f(:,[1 end]) = f(:,ibx);
    else
        f(:,[1 2 end-1 end]) = fq(:,[1 2 end-1 end]);
    end
    res_f([1 2 end-1 end],:) = 0;                                          
    if abs(Si) > abs(Pu)
        res_f(:,[1 end]) = res_f(:,ibx);
    else
        res_f(:,[1 2 end-1 end]) = 0;
    end
    
    
    % update temperature
    V_GrdT = flxdiv(T,U+u+UBG,W+w+WBG,h,'FROMM','adv');                    % advection
    
    Lpl_T  = (diff(T(:,ic),2,1)./h^2 + diff(T(ic,:),2,2)./h^2);            % diffusion
    
    RT(ic,ic) = -V_GrdT + Lpl_T/PeT - RctR_f(ic,ic)/St;                    % total rate of change
    
    res_T = (T-To)./dt - (theta.*RT + (1-theta).*RTo);                     % residual temperature evolution equation
    
    T     = T - res_T.*dt/10;                                              % update temperature solution
    
    T([1 end],:) = To([1 end],:);                                          % apply boundary conditions
    T(:,[1 end]) = T(:,ibx);
    res_T([1 end],:) = 0;                                                  
    res_T(:,[1 end]) = 0;
    
    
    % update major element composition
    Div_fMAJV = flxdiv((1-f).*MAJs,U +UBG,W +WBG,h,'FLXDIV','flx') ...
              + flxdiv(   f .*MAJf,uf+UBG,wf+WBG,h,'FLXDIV','flx');        % advection/compaction
    
    Lpl_MAJ   = (diff(MAJ(:,ic),2,1)./h^2 + diff(MAJ(ic,:),2,2)./h^2);     % diffusion
    
    RMAJ(ic,ic) = -Div_fMAJV + Lpl_MAJ/PeC;                                % total rate of change
    
    res_MAJ = (MAJ-MAJo)./dt - (theta.*RMAJ + (1-theta).*RMAJo);           % residual composition evolution equation
    
    MAJ     = MAJ - res_MAJ.*dt/10;                                        % update composition solution
    
    MAJ     = max(1e-16,min(1-1e-16,MAJ));                                 % enforce min/max bounds
        
    MAJ([1 end],:) = MAJo([1 end],:);                                      % apply boundary conditions
    MAJ(:,[1 end]) = MAJ(:,ibx);
    res_MAJ([1 end],:) = 0;
    res_MAJ(:,[1 end]) = 0;
    
    
    % *****  TRACE ELEMENTS  **********************************************

    % update incompatible trace element composition
    Div_fTRIV = flxdiv((1-f).*TRIs,U +UBG,W +WBG,h,'FROMM','flx') ...
                     + flxdiv(   f .*TRIf,uf+UBG,wf+WBG,h,'FROMM','flx');  % advection/compaction
    
    Lpl_TRI   = (diff(TRI(:,ic),2,1)./h^2 + diff(TRI(ic,:),2,2)./h^2);     % diffusion
    
    RTRI(ic,ic) = -Div_fTRIV + Lpl_TRI/PeC;                                % total rate of change
    
    res_TRI = (TRI-TRIo)./dt - (theta.*RTRI + (1-theta).*RTRIo);           % residual composition evolution equation
    
    TRI     = TRI - res_TRI.*dt/10;                                        % update composition solution
    
    TRI([1 end],:) = TRIo([1 end],:);                                      % apply boundary conditions
    TRI(:,[1 end]) = TRI(:,ibx);
    res_TRI([1 end],:) = 0;
    res_TRI(:,[1 end]) = 0;
    
    
    % update compatible trace element composition
    Div_fTRCV = flxdiv((1-f).*TRCs,U +UBG,W +WBG,h,'FROMM','flx') ...
              + flxdiv(   f .*TRCf,uf+UBG,wf+WBG,h,'FROMM','flx');         % advection/compaction
    
    Lpl_TRC   = (diff(TRC(:,ic),2,1)./h^2 + diff(TRC(ic,:),2,2)./h^2);     % diffusion
    
    RTRC(ic,ic) = -Div_fTRCV + Lpl_TRC/PeC;                                % total rate of change
    
    res_TRC = (TRC-TRCo)./dt - (theta.*RTRC + (1-theta).*RTRCo);           % residual composition evolution equation
    
    TRC     = TRC - res_TRC.*dt/10;                                        % update composition solution
    
    TRC([1 end],:) = TRCo([1 end],:);                                      % apply boundary conditions
    TRC(:,[1 end]) = TRC(:,ibx);
    res_TRC([1 end],:) = 0;
    res_TRC(:,[1 end]) = 0;
    

    % *****  RADIOGENIC ISOTOPES  *****************************************

    % update incompatible trace element composition
    Div_fIRPV = flxdiv((1-f).*IRPs,U +UBG,W +WBG,h,'FROMM','flx') ...
              + flxdiv(   f .*IRPf,uf+UBG,wf+WBG,h,'FROMM','flx');         % advection/compaction
    
    Lpl_IRP   = (diff(IRP(:,ic),2,1)./h^2 + diff(IRP(ic,:),2,2)./h^2);     % diffusion
    
    RIRP(ic,ic) = -Div_fIRPV + Lpl_IRP/PeC - DcyR_IR(ic,ic);               % total rate of change
    
    res_IRP = (IRP-IRPo)./dt - (theta.*RIRP + (1-theta).*RIRPo);           % residual composition evolution equation
    
    IRP     = IRP - res_IRP.*dt/10;                                        % update composition solution
    
    IRP([1 end],:) = IRP(ibz,:);                                           % apply boundary conditions
    IRP(:,[1 end]) = IRP(:,ibx);
    res_IRP([1 end],:) = 0;
    res_IRP(:,[1 end]) = 0;
    
    
    % update compatible trace element composition
    Div_fIRDV = flxdiv((1-f).*IRDs,U +UBG,W +WBG,h,'FROMM','flx') ...
              + flxdiv(   f .*IRDf,uf+UBG,wf+WBG,h,'FROMM','flx');         % advection/compaction
    
    Lpl_IRD   = (diff(IRD(:,ic),2,1)./h^2 + diff(IRD(ic,:),2,2)./h^2);     % diffusion
    
    RIRD(ic,ic) = -Div_fIRDV + Lpl_IRD/PeC + DcyR_IR(ic,ic);               % total rate of change
    
    res_IRD = (IRD-IRDo)./dt - (theta.*RIRD + (1-theta).*RIRDo);           % residual composition evolution equation
    
    IRD     = IRD - res_IRD.*dt/10;                                        % update composition solution
    
    IRD([1 end],:) = IRD(ibz,:);                                           % apply boundary conditions
    IRD(:,[1 end]) = IRD(:,ibx);
    res_IRD([1 end],:) = 0;
    res_IRD(:,[1 end]) = 0;
    
    
    % *****  STABLE ISOTOPES  *********************************************

    % update solid stable isotope composition
    V_GrdISS = flxdiv(ISS,U+UBG,W+WBG,h,'FROMM','adv');                    % advection/compaction
    
    Lpl_ISS  = (diff(ISS(:,ic),2,1)./h^2 + diff(ISS(ic,:),2,2)./h^2);      % diffusion
    
    RctR_ISS = -(ISR-ISS).*RctR_f./max(flim,1-f);                         % reactive transfer rate

    RISS(ic,ic) = -V_GrdISS + Lpl_ISS/PeC + RctR_ISS(ic,ic);               % total rate of change
    
    res_ISS = (ISS-ISSo)./dt - (theta.*RISS + (1-theta).*RISSo);           % residual composition evolution equation
    
    ISS     = ISS - res_ISS.*dt/10;                                        % update composition solution
    
    ISS([1 end],:) = ISS(ibz,:);                                           % apply boundary conditions
    ISS(:,[1 end]) = ISS(:,ibx);    
    res_ISS([1 end],:) = 0;
    res_ISS(:,[1 end]) = 0;
    
    
    % update fluid stable isotope composition
    V_GrdISF = flxdiv(ISF,uf+UBG,wf+WBG,h,'FROMM','adv');                  % advection/compaction
    
    Lpl_ISF  = (diff(ISF(:,ic),2,1)./h^2 + diff(ISF(ic,:),2,2)./h^2);      % diffusion
    
    RctR_ISF = (ISR-ISF).*RctR_f./max(flim,f);                             % reactive transfer rate
    
    RISF(ic,ic) = -V_GrdISF + Lpl_ISF/PeC + RctR_ISF(ic,ic);               % total rate of change
    
    res_ISF = (ISF-ISFo)./dt - (theta.*RISF + (1-theta).*RISFo);           % residual composition evolution equation
    
    ISF     = ISF - res_ISF.*dt/10;                                        % update composition solution
    
    ISF([1 end],:) = ISF(ibz,:);                                           % apply boundary conditions
    ISF(:,[1 end]) = ISF(:,ibx);    
    res_ISF([1 end],:) = 0;
    res_ISF(:,[1 end]) = 0;
    
end
