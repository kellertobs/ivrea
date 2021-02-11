% calculate phase equilibrium

function [f,cs,cf]  =  equilibrium(T,C,perT,perC,dperT,dperC,PhDg)

% cs = max(0,min(C,erfc((2+PhDg).*   max(0,min(1,T))) ));
% cf = max(C,min(1,erf ( 2      .*(1-max(0,min(1,T))))));

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

f = max(0,min(1,(C-cs) ./ (cf-cs))); 

end