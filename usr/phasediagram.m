clear;

T = linspace(0,1,1e3);

PhDg = 2;
perC = 0.7;
perT = 0.3;
dperC = 0.01;
dperT = 0.01;

perC = perC-dperC;
perT = perT-dperT;
cs1 = max(0,min(1, erfc((1.5+PhDg).*(T-perT)./(1-perT)).*perC ));
cs2 = max(0,min(1, perC+(1-perC).*erfc((1.5+PhDg).*(T)./perT) ));

cs = zeros(size(T));
cs(T>=perT) = cs1(T>=perT);
cs(T< perT) = cs2(T< perT);

perC = perC+2*dperC;
perT = perT+2*dperT;
cf1 = max(0,min(1, perC.*erf(2.*(1-T)./(1-perT))./erf(2) ));
cf2 = max(0,min(1, perC+(1-perC).*erf(2.*(1-T-(1-perT))./(perT))./erf(2) ));

cf = zeros(size(T));
cf(T>=perT) = cf1(T>=perT);
cf(T< perT) = cf2(T< perT);

figure(1); clf;
plot(cs1,T,'b--',cs2,T,'b-.',cf1,T,'r--',cf2,T,'r-.',cs,T,'k-','LineWidth',1.5); axis tight; hold on; box on;
plot(cf,T,'k-','LineWidth',1.5);

figure(2); clf;
plot(T,cs1,'b--',T,cs2,'b-.',T,cf1,'r--',T,cf2,'r-.',T,cs,'k-','LineWidth',1.5); axis tight; hold on; box on;
plot(T,cf,'k-','LineWidth',1.5);