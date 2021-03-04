% thermochemical solver

Div_fV(ic,ic) = flxdiv(1-f,U+UBG,W+WBG,h,ADVN,'flx');                      % phase advection/compaction

if step > 0
    
    % update melt fraction
    Lpl_f(ic,ic)  = (diff(f(:,ic),2,1)./h^2 ...
                  +  diff(f(ic,:),2,2)./h^2);                              % diffusion
    
    % RbndF  = max(0,-(f-fin)/5/dt .* exp(-(( Z)-L/2).^2./(L/20)^2)) ...
    % + min(0, (f-fin)/5/dt .* exp(-((-Z)-L/2).^2./(L/20)^2));
    
    res_f = (f-fo)./dt - (  theta .*(Div_fV  + Lpl_f /PeC + RctR_f ) ...
                       + (1-theta).*(Div_fVo + Lpl_fo/PeC + RctR_fo));     % residual phase evolution equation
    
    res_f([1 2 end-1 end],:) = 0;
    if abs(Si) > abs(Pu)
        res_f(:,[1 end]) = res_f(:,ibx);
    else
        res_f(:,[1 2 end-1 end]) = 0;
    end
    
    f = f - res_f.*dt/10;                                                  % update phase fraction
    
    f = max(1e-16,min(1-1e-16,f));                                         % enforce min/max bounds
    f(fq<=10*flim) = fq(fq<=10*flim);

    f([1 2 end-1 end],:) = fq([1 2 end-1 end],:);
    if abs(Si) > abs(Pu)
        f(:,[1 end]) = f(:,ibx);
    else
        f(:,[1 2 end-1 end]) = fq(:,[1 2 end-1 end]);
    end                                             

    % update temperature
    V_GrdT(ic,ic) = flxdiv(T,U+u+UBG,W+w+WBG,h,ADVN,'adv');                % advection
    
    Lpl_T(ic,ic)  = (diff(T(:,ic),2,1)./h^2 ...
                  +  diff(T(ic,:),2,2)./h^2);                              % diffusion
    
    % RbndT  = max(0,-(T-Tin)/5/dt .* exp(-(( Z)-L/2).^2./(L/20)^2)) ...
    % + min(0,-(T-Tin)/5/dt .* exp(-((-Z)-L/2).^2./(L/20)^2));
    
    res_T = (T-To)./dt - (  theta .*(-V_GrdT  + Lpl_T /PeT - RctR_f /St) ...
                       + (1-theta).*(-V_GrdTo + Lpl_To/PeT - RctR_fo/St)); % residual temperature evolution equation
        
    res_T([1 end],:) = 0;
    res_T(:,[1 end]) = 0;
    
    T = T - res_T.*dt/10;                                                  % update temperature solution
    
    T([1 end],:) = To([1 end],:);                                          % apply boundary conditions
    T(:,[1 end]) = T(:,ibx);
    
    % update major element composition
    Div_fMAJV(ic,ic) = flxdiv((1-f).*MAJs,U +UBG,W +WBG,h,ADVN,'flx') ...
                     + flxdiv(   f .*MAJf,uf+UBG,wf+WBG,h,ADVN,'flx');     % advection/compaction
    
    Lpl_MAJ(ic,ic)  = (diff(MAJ(:,ic),2,1)./h^2 ...
                    +  diff(MAJ(ic,:),2,2)./h^2);                          % diffusion
    
    % RbndMAJ = max(0,-(MAJ-MAJin)/5/dt .* exp(-(( Z)-L/2).^2./(L/20)^2)) ...
    % + min(0,-(MAJ-MAJin)/5/dt .* exp(-((-Z)-L/2).^2./(L/20)^2));
    
    res_MAJ = (MAJ-MAJo)./dt - (  theta .*(-Div_fMAJV  + Lpl_MAJ /PeC) ...
                             + (1-theta).*(-Div_fMAJVo + Lpl_MAJo/PeC));   % residual composition evolution equation
        
    res_MAJ([1 end],:) = 0;
    res_MAJ(:,[1 end]) = 0;
    
    MAJ = MAJ - res_MAJ.*dt/10;                                            % update composition solution
    
    MAJ = max(1e-16,min(1-1e-16,MAJ));                                     % enforce min/max bounds
        
    MAJ([1 end],:) = MAJo([1 end],:);                                      % apply boundary conditions
    MAJ(:,[1 end]) = MAJ(:,ibx);
    
    % update incompatible trace element composition
    Div_fTRIV(ic,ic) = flxdiv((1-f).*TRIs,U +UBG,W +WBG,h,ADVN,'flx') ...
                     + flxdiv(   f .*TRIf,uf+UBG,wf+WBG,h,ADVN,'flx');     % advection/compaction
    
    Lpl_TRI(ic,ic)  = (diff(TRI(:,ic),2,1)./h^2 ...
                    +  diff(TRI(ic,:),2,2)./h^2);                          % diffusion
    
    % RbndTRI = max(0,-(TRI-TRIin)/5/dt .* exp(-(( Z)-L/2).^2./(L/20)^2)) ...
    % + min(0,-(TRI-TRIin)/5/dt .* exp(-((-Z)-L/2).^2./(L/20)^2));
    
    res_TRI = (TRI-TRIo)./dt - (  theta .*(-Div_fTRIV  + Lpl_TRI /PeC) ...
                             + (1-theta).*(-Div_fTRIVo + Lpl_TRIo/PeC));   % residual composition evolution equation
    
    res_TRI([1 end],:) = 0;
    res_TRI(:,[1 end]) = 0;
    
    TRI = TRI - res_TRI.*dt/10;                                            % update composition solution
    
    TRI([1 end],:) = TRIo([1 end],:);                                      % apply boundary conditions
    TRI(:,[1 end]) = TRI(:,ibx);
    
    % update compatible trace element composition
    Div_fTRCV(ic,ic) = flxdiv((1-f).*TRCs,U +UBG,W +WBG,h,ADVN,'flx') ...
                     + flxdiv(   f .*TRCf,uf+UBG,wf+WBG,h,ADVN,'flx');     % advection/compaction
    
    Lpl_TRC(ic,ic)  = (diff(TRC(:,ic),2,1)./h^2 ...
                    +  diff(TRC(ic,:),2,2)./h^2);                          % diffusion
    
    % RbndTRC = max(0,-(TRC-TRCin)/5/dt .* exp(-(( Z)-L/2).^2./(L/20)^2)) ...
    % + min(0,-(TRC-TRCin)/5/dt .* exp(-((-Z)-L/2).^2./(L/20)^2));
    
    res_TRC = (TRC-TRCo)./dt - (  theta .*(-Div_fTRCV  + Lpl_TRC /PeC) ...
                             + (1-theta).*(-Div_fTRCVo + Lpl_TRCo/PeC));   % residual composition evolution equation
    
    res_TRC([1 end],:) = 0;
    res_TRC(:,[1 end]) = 0;
    
    TRC = TRC - res_TRC.*dt/10;                                            % update composition solution
    
    TRC([1 end],:) = TRCo([1 end],:);                                      % apply boundary conditions
    TRC(:,[1 end]) = TRC(:,ibx);
    
    % update solid stable isotope composition
    Div_fISSV(ic,ic) = flxdiv(ISS,U+UBG,W+WBG,h,ADVN,'adv');               % advection/compaction
    
    Lpl_ISS(ic,ic)  = (diff(ISS(:,ic),2,1)./h^2 ...
                    +  diff(ISS(ic,:),2,2)./h^2);                          % diffusion
    
    % RbndISS = max(0,-(ISS-ISSin)/5/dt .* exp(-(( Z)-L/2).^2./(L/20)^2)) ...
    % + min(0,-(ISS-ISSin)/5/dt .* exp(-((-Z)-L/2).^2./(L/20)^2));
    
    res_ISS = (ISS-ISSo)./dt - (  theta .*(-Div_fISSV  + Lpl_ISS /PeC - (ISR -ISS ).*RctR_f ./max(flim,1-f )) ...
                             + (1-theta).*(-Div_fISSVo + Lpl_ISSo/PeC - (ISRo-ISSo).*RctR_fo./max(flim,1-fo)));  % residual composition evolution equation
    
    res_ISS([1 end],:) = 0;
    res_ISS(:,[1 end]) = 0;
    
    ISS = ISS - res_ISS.*dt/10;                                            % update composition solution
    
    ISS([1 end],:) = ISSo([1 end],:);                                      % apply boundary conditions
    ISS(:,[1 end]) = ISS(:,ibx);

    % update fluid stable isotope composition
    Div_fISFV(ic,ic) = flxdiv(ISF,uf+UBG,wf+WBG,h,ADVN,'adv');             % advection/compaction
    
    Lpl_ISF(ic,ic)  = f(ic,ic).* (diff(ISF(:,ic),2,1)./h^2 ...
                               +  diff(ISF(ic,:),2,2)./h^2);               % diffusion
    
    % RbndISF = max(0,-(ISF-ISFin)/5/dt .* exp(-(( Z)-L/2).^2./(L/20)^2)) ...
    % + min(0,-(ISF-ISFin)/5/dt .* exp(-((-Z)-L/2).^2./(L/20)^2));
    
    res_ISF = (ISF-ISFo)./dt - (  theta .*(-Div_fISFV  + Lpl_ISF /PeC + (ISR -ISF ).*RctR_f ./max(flim,f )) ...
                             + (1-theta).*(-Div_fISFVo + Lpl_ISFo/PeC + (ISRo-ISFo).*RctR_fo./max(flim,fo)));  % residual composition evolution equation
    
    res_ISF([1 end],:) = 0;
    res_ISF(:,[1 end]) = 0;
    
    ISF = ISF - res_ISF.*dt/10;                                            % update composition solution
    
    ISF([1 end],:) = ISFo([1 end],:);                                      % apply boundary conditions
    ISF(:,[1 end]) = ISF(:,ibx);
    
end
