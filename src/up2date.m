%*****  TWO-PHASE RHEOLOGY  ***********************************************

% update yield stress
twophsi = double(f>=flim);
for k  = 1:5                                                   % regularisation
    kk = 1;
    twophsi(ic,ic) = twophsi(ic,ic) + kk.*(diff(twophsi(ic,:),2,2)+diff(twophsi(:,ic),2,1))./8;
    twophsi([1 end],:) = twophsi(ibz,:);
    twophsi(:,[1 end]) = twophsi(:,ibx);
end
twophs = (1-gamma).*min(1,2*twophsi) + gamma*twophs;

Pefct = (1+erf(2*(log10(f)+3)))/2;
Pe = Pefct.*p + (1-Pefct).*Pt;
yieldt_GM = max(1e-6, 1*Ty + Pe  ) .* YDMG.^DMG + etamin.*eps + bnchmrk*10;
yieldt_MC = max(1e-6, 2*Ty + Pe/2) .* YDMG.^DMG + etamin.*eps + bnchmrk*10;
yieldt    = min(yieldt_GM,yieldt_MC);% .* YDMG.^DMG + etamin.*eps + bnchmrk*10;

yieldp = min(-1e-6,((tau - etamin.*eps - bnchmrk*10) .* YDMG.^-DMG - 1*Ty - (1-Pefct).*Pt)./max(1e-3,Pefct));   
% yieldp_GM = min(-1e-6, -1*Ty +   tau );
% yieldp_MC = min(-1e-2, -4*Ty + 2*tau );

% yieldt_GM = max(1e-3, (1*Ty + Pe   + tau)/2 );
% yieldt_MC = max(1e-3, (2*Ty + Pe/2 + tau)/2);
% yieldt_SP = max(1e-3,  2*Ty + Pt/2);
% yieldt    = twophs.*min(yieldt_GM,yieldt_MC) + (1-twophs).*yieldt_MC;% + (1-twophs).*yieldt_MC);% .* YDMG.^DMG + etamin.*eps + bnchmrk*10;
% 
% yieldp_GM = max(-1*Ty+1e-3, (-1*Ty + Pe + 1*tau)/2);
% yieldp_MC = max(-4*Ty+1e-3, (-4*Ty + Pe + 2*tau)/2);
% yieldp_SP = Pt;
% yieldp    = twophs.*max(yieldp_GM,yieldp_MC) + (1-twophs).*yieldp_MC;% + (1-twophs).*yieldp_MC);% .* YDMG.^DMG;

% update viscosities
etav  = exp(Es*(1./T-1./T0_eta) - lambda.*(f-f0_eta)) .* (1/2+1/2*(max(1e-16,eps-epsELA)./e0_eta).^-n) .* EMAJ.^MAJ; 
etav  = (1./etamax + 1./etav).^-1 + etamin;
etav([1 end],:) = etav(ibz,:);
etav(:,[1 end]) = etav(:,ibx);
    
etay  =  yieldt./max(1e-16,eps-epsELA);                                                      % shear visco-plasticity
etay  =  min(etav,etay);

zetav = etay./(f .* (1-f).^0.5);
zetav([1 end],:) = zetav(ibz,:);
zetav(:,[1 end]) = zetav(:,ibx);

for k  = 1:ceil(kappa)                                                     % regularisation
    kk = kappa/ceil(kappa);
    etay(ic,ic) = etay(ic,ic) + kk.*(diff(etay(ic,:),2,2)+diff(etay(:,ic),2,1))./8;
    etay([1 end],:) = etay(ibz,:);
    etay(:,[1 end]) = etay(:,ibx);
end

zetay = -yieldp./max(1e-16,ups-upsELA);                                           % compaction visco-plasticity
zetay = min(zetav,zetay);

for k  = 1:ceil(kappa)                                                     % regularisation
    kk = kappa/ceil(kappa);
    zetay(ic,ic) = zetay(ic,ic) + kk.*(diff(zetay(ic,:),2,2)+diff(zetay(:,ic),2,1))./8;
    zetay([1 end],:) = zetay(ibz,:);
    zetay(:,[1 end]) = zetay(:,ibx);
end

 eta  =   etay.*(1-gamma) +  eta.*gamma;                                   % effective shear viscosity
zeta  =  zetay.*(1-gamma) + zeta.*gamma;                                   % effective shear viscosity
zeta  =  (flim + 1./zeta ).^-1;
zetav =  (flim + 1./zetav).^-1;

etac  = (eta(im,im)+eta(ip,im)+eta(im,ip)+eta(ip,ip)).*0.25;  % interpolate to cell corners

Gdt = De*dt + etamin/10;
Kdt = De*dt./(max(flim,f).^0.5) + etamin/10;
Kdt([1 end],:) = Kdt(ibz,:);
Kdt(:,[1 end]) = Kdt(:,ibx);

 eta_vep = (1./ eta + 1./Gdt).^-1;
zeta_vep = (1./zeta + 1./Kdt).^-1;

chi_vep  = (1 + Gdt./ eta).^-1;
 xi_vep  = (1 + Kdt./zeta).^-1;

etac_vep = (1./ etac + 1./Gdt).^-1;
chic_vep = (1 + Gdt./etac).^-1;

% etac_vep  = (eta_vep(im,im)+eta_vep(ip,im)+eta_vep(im,ip)+eta_vep(ip,ip)).*0.25;  % interpolate to cell corners
% chic_vep  = (chi_vep(im,im)+chi_vep(ip,im)+chi_vep(im,ip)+chi_vep(ip,ip)).*0.25;  % interpolate to cell corners

% update segregation coefficient
K  = (f/f0).^3 .* (1-f).^2 .* exp(-Ef*(1./T-1./T0)) .* KDMG.^DMG;  % segregation coefficient

for k  = 1:ceil(kappa)                                                     % regularisation
    kk = kappa/ceil(kappa);
    K(ic,ic) = K(ic,ic) + kk.*(diff(K(ic,:),2,2)+diff(K(:,ic),2,1))./8;
    K([1 end],:) = K(ibz,:);
    K(:,[1 end]) = K(:,ibx);
end

K  = 1./(1./K + 1e-3);
K  = (flim/f0)^3 + K;

