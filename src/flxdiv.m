function [Div_va] = flxdiv(a,u,w,h,type)

wp = w(2:end  ,2:end-1);
wm = w(1:end-1,2:end-1);
up = u(2:end-1,2:end  );
um = u(2:end-1,1:end-1);

agh                    = zeros(size(a)+2);
agh(2:end-1,2:end-1)   = a;

agh([1 2 end-1 end],:) = agh([4 3 end-2 end-3],:);
agh(:,[1 2 end-1 end]) = agh(:,[end-2 end-3 4 3]);

acc = agh(3:end-2,3:end-2);
ajp = agh(4:end-1,3:end-2);  ajpp = agh(5:end-0,3:end-2);
ajm = agh(2:end-3,3:end-2);  ajmm = agh(1:end-4,3:end-2);
aip = agh(3:end-2,4:end-1);  aipp = agh(3:end-2,5:end-0);
aim = agh(3:end-2,2:end-3);  aimm = agh(3:end-2,1:end-4);

Div_va = +     up .*(-(aipp-aip)./h./8 + (aip + acc)./h./2 + (acc-aim )./h./8) ...
         - abs(up).*(-(aipp-aip)./h./8 + (aip - acc)./h./4 - (acc-aim )./h./8) ...
         -     um .*(-(aip -acc)./h./8 + (acc + aim)./h./2 + (aim-aimm)./h./8) ...
         + abs(um).*(-(aip -acc)./h./8 + (acc - aim)./h./4 - (aim-aimm)./h./8) ...
         +     wp .*(-(ajpp-ajp)./h./8 + (ajp + acc)./h./2 + (acc-ajm )./h./8) ...
         - abs(wp).*(-(ajpp-ajp)./h./8 + (ajp - acc)./h./4 - (acc-ajm )./h./8) ...
         -     wm .*(-(ajp -acc)./h./8 + (acc + ajm)./h./2 + (ajm-ajmm)./h./8) ...
         + abs(wm).*(-(ajp -acc)./h./8 + (acc - ajm)./h./4 - (ajm-ajmm)./h./8);
                      
if strcmp(type,'adv')
    Div_va = Div_va - acc.*(diff(w(:,2:end-1),1,1)./h + diff(u(2:end-1,:),1,2)./h);
end
    
end