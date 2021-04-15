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


plimi = double(f>=flim);
for k  = 1:10                                                   % regularisation
    kk = 1;
    plimi(ic,ic) = plimi(ic,ic) + kk.*(diff(plimi(ic,:),2,2)+diff(plimi(:,ic),2,1))./8;
    plimi([1 end],:) = plimi(ibz,:);
    plimi(:,[1 end]) = plimi(:,ibx);
end
plim  = min(1,2.*plimi);
% plim  = gamma.*plim + (1-gamma).*plimi;

Pe = p.*plim + Pt.*(1-plim);
yieldt_GM = max(1e-16, 1 + Pe  );
yieldt_MC = max(1e-16, 2 + Pe/2);
yieldt    = min(yieldt_GM,yieldt_MC) .* YDMG.^DMG + etamin.*eps + bnchmrk*10;

% update rheological parameters
etav  = exp(Es*(1./T-1./T0) -lambda.*(f-f0)) .* (1/2+1/2*(eps./eps0).^-n); 
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

 eta  =   etay.*(1-gamma) +  eta.*gamma;                                   % effective shear viscosity
zeta  =  zetay.*(1-gamma) + zeta.*gamma;                                   % effective shear viscosity

etac  = (eta(im,im)+eta(ip,im) + eta(im,ip)+eta(ip,ip)).*0.25;             % interpolate to cell corners

limf  = plim.*limf + (1-plim).*f;
K     = log10((limf/f0).^3 .* (1-limf).^1.5 .* exp(-Ef*(1./T-1./T0))) + KDMG.*DMG;  % segregation coefficient

% update iterative and physical time step sizes
dtW = (10.^(( eta(im,:)+ eta(ip,:)).*0.5)./(h/2)^2 ...
    +  10.^((zeta(im,:)+zeta(ip,:)).*0.5)./(h/2)^2).^-1;                   % W iterative step size
dtU = (10.^(( eta(:,im)+ eta(:,ip)).*0.5)./(h/2)^2 ...
    +  10.^((zeta(:,im)+zeta(:,ip)).*0.5)./(h/2)^2).^-1;                   % U iterative step size
dtP = (1./(10.^eta) + 10.^K./(h/2)^2).^-1;                                 % P iterative step size
