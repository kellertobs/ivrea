% update tensor magnitudes
eps(ic,ic) = (  (exx(ic,ic).^2 + ezz(ic,ic).^2 ...
           + 2.*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2).*0.25)./2).^0.5 + 1e-16;
eps([1 end],:) = eps(ibz,:);                                               % periodic boundaries
eps(:,[1 end]) = eps(:,ibx);

tau(ic,ic) = (  (txx(ic,ic).^2 + tzz(ic,ic).^2 ...
           + 2.*(txz(1:end-1,1:end-1).^2 + txz(2:end,1:end-1).^2 + txz(1:end-1,2:end).^2 + txz(2:end,2:end).^2).*0.25)./2).^0.5 + 1e-16;
tau([1 end],:) = tau(ibz,:);                                               % periodic boundaries
tau(:,[1 end]) = tau(:,ibx);

yieldt = max(1e-16,1 + p) + etamin.*eps + bnchmrk*10;

% update rheological parameters
etav  = exp(Es*(1./T-1./T0)) .* (1-(f-f0)./0.3).^8 .* (1/2+1/2*(eps./eps0).^((1-n)./n));
etav  = (1./etamax + 1./etav).^-1 + etamin;

etay  =  log10(yieldt)-log10(eps);                                         % shear visco-plasticity
etay  =  min(log10(etav),etay);

for k  = 1:ceil(kappa)                                                     % regularisation
    kk = kappa/ceil(kappa);
    etay(ic,ic) = etay(ic,ic) + kk.*(diff(etay(ic,:),2,2)+diff(etay(:,ic),2,1))./8;
    etay([1 end],:) = etay(ibz,:);
    etay(:,[1 end]) = etay(:,ibx);
end

eta  =  etay.*(1-delta) + eta.*delta;                                      % effective shear viscosity

etac = (eta(im,im)+eta(ip,im) + eta(im,ip)+eta(ip,ip)).*0.25;              % interpolate to cell corners

zeta = eta - log10(max(1e-2,f.*(1-f).^0.5));                     % cmpct viscosity
zeta = max(log10(etamin/f0),zeta);

K    = (f/f0).^3 .* (1-f).^1.5 .* exp(-Ef*(1./T-1./T0));              % segregation coefficient
            
% update iterative and physical time step sizes
dtW = (10.^(( eta(im,:)+ eta(ip,:)).*0.5)./(h/2)^2 ...
    +  10.^((zeta(im,:)+zeta(ip,:)).*0.5)./(h/2)^2).^-1;                   % W iterative step size
dtU = (10.^(( eta(:,im)+ eta(:,ip)).*0.5)./(h/2)^2 ...
    +  10.^((zeta(:,im)+zeta(:,ip)).*0.5)./(h/2)^2).^-1;                   % U iterative step size
dtP = 0.25.*(1./(10.^eta) + K./(h/2)^2).^-1;                                     % P iterative step size

if step > 0
    % physical time step
    Vel = [U(:)+UBG(:);W(:)+WBG(:);u(:);w(:);uf(:);wf(:)];                     % combine all velocity components
    dt  = CFL*min([h/2/max(abs(Vel)), 0.005./max(abs(Div_fV(:)))]);            % physical time step

    % update phase equilibrium
    [fq,MAJsq,MAJfq]  =  equilibrium(T,MAJ,perT,perCs,perCf,PhDg);
    RctR_f            =  0.75.*RctR_f + 0.25.*Da.*(fq-f);
    for k  = 1:ceil(kappa)                                                     % regularisation
        kk = kappa/ceil(kappa);
        RctR_f(ic,ic) = RctR_f(ic,ic) + kk.*(diff(RctR_f(ic,:),2,2)+diff(RctR_f(:,ic),2,1))./8;
        RctR_f([1 end],:) = 0;
        RctR_f(:,[1 end]) = RctR_f(:,ibx);
    end
    
    % update phase major element composition
    KMAJ = MAJsq./MAJfq;
    MAJf = MAJ./(f + (1-f).*KMAJ);
    MAJs = MAJ./(f./KMAJ + (1-f));
    MAJf(f<=  1e-3) = MAJfq(f<=1e-3);  MAJs(f<=  1e-3) = MAJ(f<=    1e-3);
    MAJf(f>=1-1e-3) = MAJ(f>=1-1e-3);  MAJs(f>=1-1e-3) = MAJsq(f>=1-1e-3);
    
    % update phase trace element composition
    TRCf = TRC./(f + (1-f).*KTRC);
    TRCs = TRC./(f./KTRC + (1-f));
    TRCs(f<=  1e-3) = TRC(f<=  1e-3);
    TRCf(f>=1-1e-3) = TRC(f>=1-1e-3);
end