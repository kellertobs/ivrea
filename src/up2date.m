%*****  TWO-PHASE RHEOLOGY  ***********************************************

% update yield stress
twophs = double(f>=flim);
for k  = 1:5                                                   % regularisation
    kk = 1;
    twophs(ic,ic) = twophs(ic,ic) + kk.*(diff(twophs(ic,:),2,2)+diff(twophs(:,ic),2,1))./8;
    twophs([1 end],:) = twophs(ibz,:);
    twophs(:,[1 end]) = twophs(:,ibx);
end
twophs  = min(1,2*twophs);

Pe = p.*twophs + Pt.*(1-twophs);
yieldt_GM = max(1e-4, 1*Ty + Pe  );
yieldt_MC = max(1e-4, 2*Ty + Pe/2);
yieldt    = (twophs.*min(yieldt_GM,yieldt_MC) + (1-twophs).*yieldt_MC) .* YDMG.^DMG + etamin.*eps + bnchmrk*10;

yieldp    = min(-1e-4, -1*Ty + tau  );
yieldp    = yieldp .* YDMG.^DMG;

% update viscosities
etav  = exp(Es*(1./T-1./T0_eta) - lambda.*(f-f0_eta)) .* (1/2+1/2*(eps./e0_eta).^-n) .* EMAJ.^MAJ; 
etav  = (1./etamax + 1./etav).^-1 + etamin;
    
etay  =  yieldt./max(1e-16,eps-epsELA);                                                      % shear visco-plasticity
etay  =  min(etav,etay);

zetav = etav./(f .* (1-f).^0.5);
zetay = -yieldp./max(1e-16,ups-upsELA);                                           % compaction visco-plasticity
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
zeta  =  min(1/flim,zeta);

etac  = (eta(im,im)+eta(ip,im)+eta(im,ip)+eta(ip,ip)).*0.25;  % interpolate to cell corners

 eta_vep = (1./ eta + 1./(De*dt+etamin/10)).^-1;
zeta_vep = (1./zeta + 1./(De*dt+etamin/10)).^-1;

chi_vep  = (1 + (De*dt+etamin/10)./ eta).^-1;
 xi_vep  = (1 + (De*dt+etamin/10)./zeta).^-1;

eta_vepc = (1./ etac + 1./(De*dt+etamin/10)).^-1;
chi_vepc = (1 + (De*dt+etamin/10)./ etac).^-1;
 
% eta_vepc  = (eta_vep(im,im)+eta_vep(ip,im)+eta_vep(im,ip)+eta_vep(ip,ip)).*0.25;  % interpolate to cell corners
% chi_vepc  = (chi_vep(im,im)+chi_vep(ip,im)+chi_vep(im,ip)+chi_vep(ip,ip)).*0.25;  % interpolate to cell corners

% update segregation coefficient
K  = (f/f0).^3 .* (1-f).^2 .* exp(-Ef*(1./T-1./T0)) .* KDMG.^DMG;  % segregation coefficient
K  = 1./(1./K + 1e-3);

for k  = 1:ceil(kappa)                                                     % regularisation
    kk = kappa/ceil(kappa);
    K(ic,ic) = K(ic,ic) + kk.*(diff(K(ic,:),2,2)+diff(K(:,ic),2,1))./8;
    K([1 end],:) = K(ibz,:);
    K(:,[1 end]) = K(:,ibx);
end

K    = max((flim/f0)^3,K);

