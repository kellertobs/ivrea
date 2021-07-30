%*****  TWO-PHASE RHEOLOGY  ***********************************************

% update yield stress
twophsi = double(f>=flim);
for k  = 1:5  % regularisation
    kk = 1;
    twophsi(ic,ic) = twophsi(ic,ic) + kk.*(diff(twophsi(ic,:),2,2)+diff(twophsi(:,ic),2,1))./8;
    twophsi([1 end],:) = twophsi(ibz,:);
    twophsi(:,[1 end]) = twophsi(:,ibx);
end
if step == 0 && it == 0
    twophs = min(1,2*twophsi);
else
    twophs = (1-gamma).*min(1,2*twophsi) + gamma*twophs;
end

Pefct     = (1+erf(2*(log10(f)+3)))/2;
Pe        = Pefct.*p + (1-Pefct).*Pt;
yieldt_GM = max(1e-6, 1*Ty + Pe  ) .* YDMG.^DMG + etamin.*eps + bnchmrk*10;
yieldt_MC = max(1e-6, 2*Ty + Pe/2) .* YDMG.^DMG + etamin.*eps + bnchmrk*10;
yieldt    = min(yieldt_GM,yieldt_MC);

yieldp    = min(-1e-6,((tau - etamin.*eps - bnchmrk*10) .* YDMG.^-DMG - 1*Ty - (1-Pefct).*Pt)./max(1e-3,Pefct));   

% update shear viscosity
etav  = exp(Es*(1./T-1./T0_eta) - lambda.*(f-f0_eta)) .* (1/2+1/2*(max(1e-16,eps-epsELA)./e0_eta).^n).^-1 .* EMAJ.^(MAJ-MAJ0_eta); 
etav  = (1./etamax + 1./etav).^-1 + etamin;
etav([1 end],:) = etav(ibz,:);
etav(:,[1 end]) = etav(:,ibx);

% update compaction viscosity
zetav = etav./(f .* (1-f).^0.5);
zetav = (flim + 1./zetav).^-1;
zetav([1 end],:) = zetav(ibz,:);
zetav(:,[1 end]) = zetav(:,ibx);

% update shear and compaction elastic moduli
Gdt = De*dt + etamin/10;
Kdt = De*dt./(max(flim,f).^0.5) + etamin/f0/10;
Kdt([1 end],:) = Kdt(ibz,:);
Kdt(:,[1 end]) = Kdt(:,ibx);
 
% update shear and compaction visco-plasticities
 etay  =   yieldt ./ max(1e-6,eps - yieldt./ etav - (yieldt-tauo)./Gdt) + etamin   ;
zetay  = - yieldp ./ max(1e-6,ups + yieldp./zetav + (yieldp-po  )./Kdt) + etamin/f0;

for k  = 1:ceil(kappa)  % regularisation
    kk = kappa/ceil(kappa);
    etay(ic,ic) = 10.^(log10(etay(ic,ic)) + kk.*(diff(log10(etay(ic,:)),2,2)+diff(log10(etay(:,ic)),2,1))./8);
    etay([1 end],:) = etay(ibz,:);
    etay(:,[1 end]) = etay(:,ibx);
end

for k  = 1:ceil(kappa)  % regularisation
    kk = kappa/ceil(kappa);
    zetay(ic,ic) = 10.^(log10(zetay(ic,ic)) + kk.*(diff(log10(zetay(ic,:)),2,2)+diff(log10(zetay(:,ic)),2,1))./8);
    zetay([1 end],:) = zetay(ibz,:);
    zetay(:,[1 end]) = zetay(:,ibx);
end

% iterative lagging
if step == 0 && it == 0
     eta_vep = (1./ etav + 1./Gdt + 1./ etay).^-1;
    zeta_vep = (1./zetav + 1./Kdt + 1./zetay).^-1;
else
     eta_vep =  eta_vep.*gamma + (1-gamma).* (1./ etav + 1./Gdt + 1./ etay).^-1;
    zeta_vep = zeta_vep.*gamma + (1-gamma).* (1./zetav + 1./Kdt + 1./zetay).^-1;
end

% visco-elastic parameters
chi_vep  =  eta_vep./Gdt;
 xi_vep  = zeta_vep./Kdt;
 
 eta =  eta_vep;
zeta = zeta_vep;

% interpolate to cell corners
etavc     = (etav   (im,im)+etav   (ip,im)+etav   (im,ip)+etav   (ip,ip)).*0.25;
etayc     = (etay   (im,im)+etay   (ip,im)+etay   (im,ip)+etay   (ip,ip)).*0.25;
etac_vep  = (eta_vep(im,im)+eta_vep(ip,im)+eta_vep(im,ip)+eta_vep(ip,ip)).*0.25;
chic_vep  = (chi_vep(im,im)+chi_vep(ip,im)+chi_vep(im,ip)+chi_vep(ip,ip)).*0.25;

% update segregation coefficient
K  = (f/f0).^3 .* (1-f).^2 .* exp(-Ef*(1./T-1./T0)) .* KDMG.^DMG;
K  = 1./(1./K + 1e-3);
K  = (flim/f0)^3 + K;
K([1 end],:) = K(ibz,:);
K(:,[1 end]) = K(:,ibx);

