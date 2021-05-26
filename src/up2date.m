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

% update yield stress
twophs = double(f>=flim);
for k  = 1:5                                                   % regularisation
    kk = 1;
    twophs(ic,ic) = twophs(ic,ic) + kk.*(diff(twophs(ic,:),2,2)+diff(twophs(:,ic),2,1))./8;
    twophs([1 end],:) = twophs(ibz,:);
    twophs(:,[1 end]) = twophs(:,ibx);
end
twophs  = min(1,2.*twophs);

Pe = p.*twophs + Pt.*(1-twophs);
yieldt_GM = max(1e-16, 1*Ty + Pe  );
yieldt_MC = max(1e-16, 2*Ty + Pe/2);
yieldt    = min(yieldt_GM,yieldt_MC) .* YDMG.^DMG + etamin.*eps + bnchmrk*10;

yieldp    = min(-1e-16, -1*Ty + tau  );
yieldp    = yieldp .* YDMG.^DMG;

% update viscosities
etav  = exp(Es*(1./T-1./T0) - lambda.*f) .* (1/2+1/2*(eps./eps0).^-n) .* EMAJ.^MAJ; 
etav  = (1./etamax + 1./etav).^-1 + etamin;
etav  = log10(etav);
etav([1 end],:) = etav(ibz,:);
etav(:,[1 end]) = etav(:,ibx);
    
etay  =  log10(yieldt)-log10(eps);                                         % shear visco-plasticity
etay  =  min(etav,etay);

% zetay =  etay - log10(f .* (1-f).^0.5);
zetav = etay - log10(f .* (1-f).^0.5);
zetay = log10(-yieldp)-log10(max(1e-16,ups));                              % compaction visco-plasticity
zetay = twophs.*min(zetav,zetay) + (1-twophs).*zetav;

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

 eta  =   etay.*(1-gamma) +  eta.*gamma;                                   % effective shear viscosity
zeta  =  zetay.*(1-gamma) + zeta.*gamma;                                   % effective shear viscosity

etac  = (eta(im,im)+eta(ip,im) + eta(im,ip)+eta(ip,ip)).*0.25;             % interpolate to cell corners

% update segregation coefficient
K  = (f/f0).^3 .* (1-f).^2 .* exp(-Ef*(1./T-1./T0)) .* KDMG.^DMG;  % segregation coefficient
K  = 1./(1./K + 1e-3);
K  = log10(K);

for k  = 1:ceil(kappa)                                                     % regularisation
    kk = kappa/ceil(kappa);
    K(ic,ic) = K(ic,ic) + kk.*(diff(K(ic,:),2,2)+diff(K(:,ic),2,1))./8;
    K([1 end],:) = K(ibz,:);
    K(:,[1 end]) = K(:,ibx);
end

zeta = min(log10(1/flim),zeta);
K    = max(log10((flim/f0)^3),K);

% update iterative and physical time step sizes
dtW = (10.^(( eta(im,:)+ eta(ip,:)).*0.5)./(h/2)^2 ...
    +  10.^((zeta(im,:)+zeta(ip,:)).*0.5)./(h/2)^2).^-1;                   % W iterative step size
dtU = (10.^(( eta(:,im)+ eta(:,ip)).*0.5)./(h/2)^2 ...
    +  10.^((zeta(:,im)+zeta(:,ip)).*0.5)./(h/2)^2).^-1;                   % U iterative step size
dtP = (1./(10.^eta) + 10.^K./(h/2)^2).^-1;                                 % P iterative step size
