function [Div_va] = flxdiv(a,u,w,h,X,Z,dr,type)

wp = w(2:end  ,2:end-1);
wm = w(1:end-1,2:end-1);
wc = (w(2:end  ,:)+w(1:end-1,:))./2;
up = u(2:end-1,2:end  );
um = u(2:end-1,1:end-1);
uc = (u(:,2:end  )+u(:,1:end-1))./2;

agh                    = zeros(size(a)+2);
agh(2:end-1,2:end-1)   = a;

% agh([1 2 end-1 end],:) = agh([end-3 end-2 3 4],:);
% agh(:,[1 2 end-1 end]) = agh(:,[end-3 end-2 3 4]);

% agh([1 2 end-1 end],:) = agh([4 3 end-2 end-3],:);
% agh(:,[1 2 end-1 end]) = agh(:,[4 3 end-2 end-3]);
% 
agh([1 2 end-1 end],:) = mean(a(:));
agh(:,[1 2 end-1 end]) = mean(a(:));

% N = length(a);
% 
% Xbnd = round(X(:,end-1)./h)    ;  Xbnd(Xbnd<2) = Xbnd(Xbnd<2)+(N-3);  Xbnd(Xbnd>N-2) = Xbnd(Xbnd>N-2)-(N-3);
% Zbnd = (1:N).';%round(Z(:,end-1)./h)+N/2;  Zbnd(Zbnd<2) = Zbnd(Zbnd<2)+(N-3);  Zbnd(Zbnd>N-2) = Zbnd(Zbnd>N-2)-(N-3);
% 
% % left boundary inflow
% idx = sub2ind(size(dr), repmat(Zbnd,1,2), Xbnd+[1 2]);
% fct = repmat(double(uc(:,  1)>0),1,2).*[1.0,0.5];
% agh(2:end-1,[    1   2]) = (1-fct).*agh(2:end-1,[    1   2]) + fct.*(mean(a(:))+dr(idx));
% 
% % right boundary inflow
% idx = sub2ind(size(dr), repmat(Zbnd,1,2), Xbnd+[-1 0]);
% fct = repmat(double(uc(:,end)<0),1,2).*[0.5,1.0];
% agh(2:end-1,[end-1 end]) = (1-fct).*agh(2:end-1,[end-1 end]) + fct.*(mean(a(:))+dr(idx));
% 
% Xbnd = (1:N);%round(X(end-1,:)./h)+N/2;  Xbnd(Xbnd<2) = Xbnd(Xbnd<2)+(N-3);  Xbnd(Xbnd>N-2) = Xbnd(Xbnd>N-2)-(N-3);
% Zbnd = round(Z(end-1,:)./h)    ;  Zbnd(Zbnd<2) = Zbnd(Zbnd<2)+(N-3);  Zbnd(Zbnd>N-2) = Zbnd(Zbnd>N-2)-(N-3);
% 
% % top boundary inflow
% idx = sub2ind(size(dr), Zbnd+[1 2].',repmat(Xbnd,2,1));
% fct = repmat(double(wc(  1,:)>0),2,1).*[1.0;0.5];
% agh([    1   2],2:end-1) = (1-fct).*agh([    1   2],2:end-1) + fct.*(mean(a(:))+dr(idx));
% 
% % bottom boundary inflow
% idx = sub2ind(size(dr), Zbnd+[-1 0].',repmat(Xbnd,2,1));
% fct = repmat(double(wc(end,:)<0),2,1).*[0.5;1.0];
% agh([end-1 end],2:end-1) = (1-fct).*agh([end-1 end],2:end-1) + fct.*(mean(a(:))+dr(idx));

acc = agh(3:end-2,3:end-2);
ajp = agh(4:end-1,3:end-2);  ajpp = agh(5:end-0,3:end-2);
ajm = agh(2:end-3,3:end-2);  ajmm = agh(1:end-4,3:end-2);
aip = agh(3:end-2,4:end-1);  aipp = agh(3:end-2,5:end-0);
aim = agh(3:end-2,2:end-3);  aimm = agh(3:end-2,1:end-4);

Div_va   =     up .*(-(aipp-aip)./h./8 + (aip + acc)./h./2 + (acc-aim )./h./8) ...
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