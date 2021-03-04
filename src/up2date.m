% update tensor magnitudes
eps(ic,ic) = (  (exx(ic,ic).^2 + ezz(ic,ic).^2 ...
           + 2.*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2).*0.25)./2).^0.5 + 1e-16;
eps([1 end],:) = eps(ibz,:);                                               % periodic boundaries
eps(:,[1 end]) = eps(:,ibx);

tau(ic,ic) = (  (txx(ic,ic).^2 + tzz(ic,ic).^2 ...
           + 2.*(txz(1:end-1,1:end-1).^2 + txz(2:end,1:end-1).^2 + txz(1:end-1,2:end).^2 + txz(2:end,2:end).^2).*0.25)./2).^0.5 + 1e-16;
tau([1 end],:) = tau(ibz,:);                                               % periodic boundaries
tau(:,[1 end]) = tau(:,ibx);

% Pe = p.*(f>=flim) + (10 + 3*B*Z + P).*(f<flim);
Pe = p;
yieldt_GM = max(1e-16, 1 + Pe  ) + etamin.*eps + bnchmrk*10; %(1-f).*0.1.^(MAJ-MAJ0)
yieldt_MC = max(1e-16, 4 + Pe/2) + etamin.*eps + bnchmrk*10; %(1-f).*0.1.^(MAJ-MAJ0)
yieldt    = yieldt_GM; %min(yieldt_GM,yieldt_MC);

% update rheological parameters
etav  = exp(Es*(1./T-1./T0) -30.*(f-f0)) .* (1/2+1/2*(eps./eps0).^((1-n)./n)); %.* 0.01.^(MAJ-MAJ0)
etav  = (1./etamax + 1./etav).^-1 + etamin;

etay  =  log10(yieldt)-log10(eps);                                         % shear visco-plasticity
etay  =  min(log10(etav),etay);

for k  = 1:ceil(kappa)                                                     % regularisation
    kk = kappa/ceil(kappa);
    etay(ic,ic) = etay(ic,ic) + kk.*(diff(etay(ic,:),2,2)+diff(etay(:,ic),2,1))./8;
    etay([1 end],:) = etay(ibz,:);
    etay(:,[1 end]) = etay(:,ibx);
end

eta  =  etay.*0.01 + eta.*0.99;                                            % effective shear viscosity

etac = (eta(im,im)+eta(ip,im) + eta(im,ip)+eta(ip,ip)).*0.25;              % interpolate to cell corners

limf = max(flim,min(1-flim,f));

zeta = eta - log10(limf.*(1-limf).^0.5);                                   % cmpct viscosity
zeta = max(log10(etamin/(10*f0)),zeta);

K    = log10((limf/f0).^3 .* (1-limf).^1.5 .* exp(-Ef*(1./T-1./T0)));      % segregation coefficient

K([1 end],:) = K(ibz,:);
K(:,[1 end]) = K(:,ibx);

% update iterative and physical time step sizes
dtW = (10.^(( eta(im,:)+ eta(ip,:)).*0.5)./(h/2)^2 ...
    +  10.^((zeta(im,:)+zeta(ip,:)).*0.5)./(h/2)^2).^-1;                   % W iterative step size
dtU = (10.^(( eta(:,im)+ eta(:,ip)).*0.5)./(h/2)^2 ...
    +  10.^((zeta(:,im)+zeta(:,ip)).*0.5)./(h/2)^2).^-1;                   % U iterative step size
dtP = 0.5*(1./(10.^eta) + 10.^K./(h/2)^2).^-1;                             % P iterative step size

if step > 0
    % physical time step
    Vel = [U(:)+UBG(:);W(:)+WBG(:);u(:);w(:);uf(:);wf(:)];           % combine all velocity components
    dt  = CFL*min([h/2/max(abs(Vel)), 0.001./max(abs(Div_fV(:)))]);        % physical time step

    % update phase equilibrium
    [fq,MAJsq,MAJfq]  =  equilibrium(T,MAJ,perT,perCs,perCf,PhDg);
    RctR_fi = min(Da,1/(2*dt)).*(fq-f);
    for k  = 1:ceil(kappa)                                                 % regularisation
        kk = kappa/ceil(kappa);
        RctR_fi(ic,ic) = RctR_fi(ic,ic) + kk.*(diff(RctR_fi(ic,:),2,2)+diff(RctR_fi(:,ic),2,1))./8;
        RctR_fi([1 end],:) = RctR_fi(ibz,:);
        RctR_fi(:,[1 end]) = RctR_fi(:,ibx);
    end
    RctR_f = 0.01.*RctR_fi + 0.99.*RctR_f;

    % update phase major element composition
    KMAJq = MAJsq./MAJfq;
    MAJf = MAJ./(f + (1-f).*KMAJq);
    MAJs = MAJ./(f./KMAJq + (1-f));
    MAJf(f<=10*flim) = MAJfq(f<=10*flim);  MAJs(f<=  flim) = MAJ(f<=    flim);
    MAJf(f>=1-flim)  = MAJ  (f>=1 -flim);  MAJs(f>=1-flim) = MAJsq(f>=1-flim);
    
    % update phase compatible trace element composition
    TRIf = TRI./(f + (1-f).*KTRI);
    TRIs = TRI./(f./KTRI + (1-f));
    TRIs(f<=  flim) = TRI(f<=  flim);
    TRIf(f>=1-flim) = TRI(f>=1-flim);
    
    % update phase compatible trace element composition
    TRCf = TRC./(f + (1-f).*KTRC);
    TRCs = TRC./(f./KTRC + (1-f));
    TRCs(f<=  flim) = TRC(f<=  flim);
    TRCf(f>=1-flim) = TRC(f>=1-flim);
    
    % stable isotope transfer composition
    ISR = double(RctR_f>=0).*ISS + double(RctR_f<0).*ISF;
end