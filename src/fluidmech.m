%% assemble coefficients for matrix velocity diagonal and right-hand side
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A
IR  = [];       % equation indeces into R
RR  = [];       % forcing entries for R

% assemble coefficients of z-stress divergence

% left boundary
ii = MapW(:,1); jj1 = ii; jj2 = MapW(:,2);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% right boundary
ii = MapW(:,end); jj1 = ii; jj2 = MapW(:,end-1);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% top boundary
ii = MapW(1,2:end-1); jj = ii;
aa = zeros(size(ii)) + 0;
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa(:)+1];
aa = -WBG(ii);
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% bottom boundary
ii = MapW(end,2:end-1); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa(:)+1];
aa = -WBG(ii);
IR = [IR; ii(:)]; RR = [RR; aa(:)];


% internal points
ii    = MapW(2:end-1,2:end-1);
EtaC1 = 10.^etac(2:end-1,1:end-1); EtaC2 = 10.^etac(2:end-1,2:end  );
EtaP1 = 10.^eta (2:end-2,2:end-1); EtaP2 = 10.^eta (3:end-1,2:end-1);

% coefficients multiplying z-velocities W
%             top          ||         bottom          ||           left            ||          right
jj1 = MapW(1:end-2,2:end-1); jj2 = MapW(3:end,2:end-1); jj3 = MapW(2:end-1,1:end-2); jj4 = MapW(2:end-1,3:end);

aa = - 2/3*(EtaP1+EtaP2)/h^2 - 1/2*(EtaC1+EtaC2)/h^2;
II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AA = [AA; aa(:)           ];      % W on stencil centre
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; 2/3*EtaP1(:)/h^2];      % W one above
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; 2/3*EtaP2(:)/h^2];      % W one below
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; 1/2*EtaC1(:)/h^2];      % W one to the left
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; 1/2*EtaC2(:)/h^2];      % W one to the right

% coefficients multiplying x-velocities U
%         top left         ||        bottom left          ||       top right       ||       bottom right
jj1 = MapU(2:end-2,1:end-1); jj2 = MapU(3:end-1,1:end-1); jj3 = MapU(2:end-2,2:end); jj4 = MapU(3:end-1,2:end);

II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; (1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % W one to the top and left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA;-(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % W one to the bottom and left
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA;-(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % W one to the top and right
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; (1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % W one to the bottom and right


% z-RHS vector
rr = zeros(size(ii)) - f0.*((f(2:end-2,2:end-1)+f(3:end-1,2:end-1))/2 - mean(f(:)));
IR = [IR; ii(:)];  RR = [RR; rr(:)];


%  assemble coefficients of x-stress divergence

% top boundary
ii = MapU(1,:); jj1 = ii; jj2 = MapU(2,:);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% bottom boundary
ii = MapU(end,:); jj1 = ii; jj2 = MapU(end-1,:);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% left side boundary
ii = MapU(2:end-1,1); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa(:)+1];
aa = -UBG(ii-NW);
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% right side boundary
ii = MapU(2:end-1,end); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa(:)+1];
aa = -UBG(ii-NW);
IR = [IR; ii(:)]; RR = [RR; aa(:)];


% internal points
ii    = MapU(2:end-1,2:end-1);
EtaC1 = 10.^etac(1:end-1,2:end-1); EtaC2 = 10.^etac(2:end  ,2:end-1);
EtaP1 = 10.^eta (2:end-1,2:end-2); EtaP2 = 10.^eta (2:end-1,3:end-1);

% coefficients multiplying x-velocities U
%            left          ||          right          ||           top             ||          bottom
jj1 = MapU(2:end-1,1:end-2); jj2 = MapU(2:end-1,3:end); jj3 = MapU(1:end-2,2:end-1); jj4 = MapU(3:end,2:end-1);

aa = - 2/3*(EtaP1+EtaP2)/h^2 - 1/2*(EtaC1+EtaC2)/h^2;
II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AA = [AA; aa(:)           ];      % U on stencil centre
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; 2/3*EtaP1(:)/h^2];      % U one to the left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; 2/3*EtaP2(:)/h^2];      % U one to the right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; 1/2*EtaC1(:)/h^2];      % U one above
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; 1/2*EtaC2(:)/h^2];      % U one below

% coefficients multiplying z-velocities W
%         top left         ||        top right          ||       bottom left       ||       bottom right
jj1 = MapW(1:end-1,2:end-2); jj2 = MapW(1:end-1,3:end-1); jj3 = MapW(2:end,2:end-2); jj4 = MapW(2:end,3:end-1);

II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; (1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % W one to the top and left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA;-(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % W one to the top and right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA;-(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % W one to the bottom and left
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; (1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % W one to the bottom and right


% x-RHS vector
rr = zeros(size(ii));
IR = [IR; ii(:)];  RR = [RR; rr(:)];


% assemble coefficient matrix & right-hand side vector
KV = sparse(II,JJ,AA,NW+NU,NW+NU);
RV = sparse(IR,ones(size(IR)),RR);

% pseudo-elastic lagging
% El = 1e6;
% KV = 1./(1./El + 1./KV);
% RV = RV - 1./(1 + El./KV).*(KV*[W(:);U(:)]);


%% assemble coefficients for segregation velocity diagonal and right-hand side
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A
IR  = [];       % equation indeces into R
RR  = [];       % forcing entries for R


% assemble coefficients of z-segregation velocity

% left boundary
ii = MapW(:,1); jj1 = ii; jj2 = MapW(:,2);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% bottom boundary
ii = MapW(:,end); jj1 = ii; jj2 = MapW(:,end-1);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% top boundary
ii = MapW(1,:).'; jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa+1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% bottom boundary
ii = MapW(end,:).'; jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa+1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];


% internal points
ii = MapW(2:end-1,2:end-1);
K1 = 10.^K(2:end-2,2:end-1); K2 = 10.^K(3:end-1,2:end-1);

% coefficients multiplying z-velocities w
II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AA = [AA; 2./(K1(:)+K2(:))];      % w on stencil centre


% z-RHS vector
rr = zeros(size(ii)) - 1;
IR = [IR; ii(:)];  RR = [RR; rr(:)];


%  assemble coefficients of x-segregation velocity

% top boundary
ii = MapU(1,:).'; jj1 = ii; jj2 = MapU(2,:).';
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% bottom boundary
ii = MapU(end,:).'; jj1 = ii; jj2 = MapU(end-1,:).';
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% left side boundary
ii = MapU(:,1); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa+1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% right side boundary
ii = MapU(:,end); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa+1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];


% internal points
ii = MapU(2:end-1,2:end-1);
K1 = 10.^K(2:end-1,2:end-2); K2 = 10.^K(2:end-1,3:end-1);

% coefficients multiplying x-velocities u
II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AA = [AA; 2./(K1(:)+K2(:))];      % u on stencil centre


% x-RHS vector
rr = zeros(size(ii));
IR = [IR; ii(:)];  RR = [RR; rr(:)];


% assemble coefficient matrix & right-hand side vector
KD = sparse(II,JJ,AA,NW+NU,NW+NU);
RD = sparse(IR,ones(size(IR)),RR);


%% assemble coefficients for gradient operator
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A


% coefficients for z-gradient
ii = MapW(2:end-1,2:end-1);

%         top              ||          bottom
jj1 = MapP(2:end-2,2:end-1); jj2 = MapP(3:end-1,2:end-1);

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)-1/h];     % one to the top
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+1/h];     % one to the bottom


% coefficients for x-gradient
ii = MapU(2:end-1,2:end-1);

%         left             ||           right
jj1 = MapP(2:end-1,2:end-2); jj2 = MapP(2:end-1,3:end-1);

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)-1/h];     % one to the left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+1/h];     % one to the right


% assemble coefficient matrix
GG = sparse(II,JJ,AA,NW+NU,NP);


%% assemble coefficients for divergence operator
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A

%internal points
ii = MapP(2:end-1,2:end-1);

% coefficients multiplying velocities U, W
%          left U          ||           right U       ||           top W           ||          bottom W
jj1 = MapU(2:end-1,1:end-1); jj2 = MapU(2:end-1,2:end); jj3 = MapW(1:end-1,2:end-1); jj4 = MapW(2:end,2:end-1);

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)-1/h];  % U one to the left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+1/h];  % U one to the right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; aa(:)-1/h];  % W one above
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; aa(:)+1/h];  % W one below

% assemble coefficient matrix
DD = sparse(II,JJ,AA,NP,NW+NU);


%% assemble coefficients for matrix pressure diagonal and right-hand side
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A
IR  = [];       % equation indeces into R
RR  = [];       % forcing entries for R


% boundary points
ii  = [MapP(1,:).'; MapP(end  ,:).']; % top & bottom
jj1 = ii;
jj2 = [MapP(2,:).'; MapP(end-1,:).'];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

ii  = [MapP(:,1); MapP(:,end  )]; % left & right
jj1 = ii;
jj2 = [MapP(:,2); MapP(:,end-1)];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];


% internal points
ii = MapP(2:end-1,2:end-1);

% coefficients multiplying matrix pressure P
% aa = (twophs(2:end-1,2:end-1).*delta/10 + (1-twophs(2:end-1,2:end-1)).*delta)./(10.^etav(2:end-1,2:end-1));
aa = delta*h^2./(10.^eta(2:end-1,2:end-1));
% aa = zeros(size(ii)) + delta;
II = [II; ii(:)]; JJ = [JJ; ii(:)];    AA = [AA; aa(:)];  % P on stencil centre


% RHS
rr = zeros(size(ii));
IR = [IR; ii(:)];
RR = [RR; rr(:)];


% assemble coefficient matrix and right-hand side vector
KP = sparse(II,JJ,AA,NP,NP);
RP = sparse(IR,ones(size(IR)),RR,NP,1);


%% assemble coefficients for compaction pressure diagonal and right-hand side
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A
IR  = [];       % equation indeces into R
RR  = [];       % forcing entries for R


% boundary points
ii  = [MapP(1,:).'; MapP(end  ,:).']; % top & bottom
jj1 = ii;
jj2 = [MapP(2,:).'; MapP(end-1,:).'];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

ii  = [MapP(:,1); MapP(:,end  )]; % left & right
jj1 = ii;
jj2 = [MapP(:,2); MapP(:,end-1)];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];


%internal points
ii = MapP(2:end-1,2:end-1);
ze = 10.^zeta(2:end-1,2:end-1);

% coefficients multiplying compaction pressure p
II = [II; ii(:)]; JJ = [JJ; ii(:)];    AA = [AA; 1./ze(:)];  % p on stencil centre

% RHS
rr = zeros(size(ii));
IR = [IR; ii(:)];
RR = [RR; rr(:)];


% assemble coefficient matrix and right-hand side vector
KC = sparse(II,JJ,AA,NP,NP);
RC = sparse(IR,ones(size(IR)),RR,NP,1);


%% assemble global coefficient matrix and right-hand side vector
OV = sparse(NU+NW,NU+NW);
OG = sparse(NU+NW,NP   );
OP = sparse(NP   ,NP   );

A = [-KV   OV    GG   GG  ; ...
      OV   KD    GG   OG  ; ...
      DD   DD    KP   OP  ; ...
      DD   OG.'  OP   KC ];

R = [RV; RD; RP; RC];



%% set w, u, p to zero where no melt
bc_ind  = [];
bc_val  = [];

twophsw = (twophs(1:end-1,:) + twophs(2:end,:))./2;
twophsu = (twophs(:,1:end-1) + twophs(:,2:end))./2;
twophsv = [twophsw(:);twophsu(:)];

ind     =  find(twophsv(:)<=1/4) + NW+NU;
bc_ind  =  [bc_ind;ind];
bc_val  =  [bc_val;zeros(size(bc_ind))];

ind     =  find(twophs (:)<=1/4) + 2*NW+2*NU+NP;
bc_ind  =  [bc_ind;ind];
bc_val  =  [bc_val;zeros(size(bc_ind))];

% ind     =  2*NW+2*NU+MapP(end-1,N/2);
% bc_ind  =  [bc_ind;ind];
% bc_val  =  [bc_val;zeros(size(bc_ind))];

% assemble and sort all boundary indices and values
[BC.ind,ind]    =  sort(bc_ind);
 BC.val         =  bc_val(ind);

% extract boundary conditions to reduce problem size
BC.free         =  1:length(A(1,:));
BC.free(BC.ind) =  [];
TMP             =  A(:,BC.ind);
R               =  R - TMP*BC.val;
A               =  A(BC.free,BC.free);
R               =  R(BC.free);


%% Scale system of equations (diagonal preconditioning)
C  =  sqrt(abs(diag(A)));
C  =  diag(sparse(1./C));

A  =  C*A*C;
R  =  C*R;


%% get residual
% get non-linear residual
F          = zeros(size(S));
F(BC.free) = A*(C\S(BC.free)) - R;

% map residual vector to 2D arrays
FW  = full(reshape(F(MapW(:))          ,(N-1), N   ));  	% z-velocity
FU  = full(reshape(F(MapU(:))          , N   ,(N-1)));  	% x-velocity
FP  = full(reshape(F(MapP(:)+2*(NW+NU)), N   , N   ));     % dynamic pressure
Fw  = full(reshape(F(MapW(:)+   NW+NU    ),(N-1), N   )); 	% z-velocity
Fu  = full(reshape(F(MapU(:)+   NW+NU    ), N   ,(N-1))); 	% x-velocity
Fp  = full(reshape(F(MapP(:)+2*(NW+NU)+NP), N   , N   )); 	% dynamic pressure
    
% report convergence
report;
        

%% Solve linear system of equations for vx, vz, P
% S(BC.free) = S(BC.free) - C*(A\F(BC.free));  % update solution
S(BC.free) = C*(A\R);  % update solution
S(BC.ind ) = BC.val;  % fill in boundary conditions

% Read out solution
% map solution vector to 2D arrays
W  = full(reshape(S(MapW(:))          ,(N-1), N   ));              % matrix z-velocity
U  = full(reshape(S(MapU(:))          , N   ,(N-1)));              % matrix x-velocity
P  = full(reshape(S(MapP(:)+2*(NW+NU)), N   , N   ));              % matrix dynamic pressure

w  = full(reshape(S(MapW(:)+   NW+NU    ),(N-1), N   ));           % segregation z-velocity
u  = full(reshape(S(MapU(:)+   NW+NU    ), N   ,(N-1)));           % segregation x-velocity
p  = full(reshape(S(MapP(:)+2*(NW+NU)+NP), N   , N   ));           % compaction pressure


%% get auxiliary variables
Wf  = W + w./max(flim,(f(im,:)+f(ip,:))./2);                       % pore fluid z-velocity
Uf  = U + u./max(flim,(f(:,im)+f(:,ip))./2);                       % pore fluid x-velocity
fWf = W.*(f(im,:)+f(ip,:))./2 + w;                                 % pore fluid z-velocity
fUf = U.*(f(:,im)+f(:,ip))./2 + u;                                 % pore fluid x-velocity
fWs = W.*(1-(f(im,:)+f(ip,:))./2);                                 % solid matrix z-velocity
fUs = U.*(1-(f(:,im)+f(:,ip))./2);                                 % solid matrix x-velocity
Pt  = Pc + P + p + 5*Z;                                            % total pressure

% update velocity divergences
ups(ic,ic) = diff(U(ic,:),1,2)./h + diff(W(:,ic),1,1)./h;          % matrix velocity divergence

ups([1 end],:) = ups(ibz,:);                                       % apply boundary conditions
ups(:,[1 end]) = ups(:,ibx);

upss(ic,ic) = diff(u(ic,:),1,2)./h + diff(w(:,ic),1,1)./h;         % segregation velocity divergence
upss([1 end],:) = upss(ibz,:);                                     % apply boundary conditions
upss(:,[1 end]) = upss(:,ibx);

% update strain rates
exx(:,ic)        = diff(U,1,2)./h - ups(:,ic)./3 - 0.*Pu;             % x-normal strain rate
exx([1 end],:)   = exx(ibz,:);                                     % apply boundary conditions
exx(:,[1 end])   = exx(:,ibx);
ezz(ic,:)        = diff(W,1,1)./h - ups(ic,:)./3 + 0.*Pu;             % z-normal strain rate
ezz([1 end],:)   = ezz(ibz,:);                                     % apply boundary conditions
ezz(:,[1 end])   = ezz(:,ibx);
exz              = 1/2.*(diff(U,1,1)./h+diff(W,1,2)./h) - 0.*Si;      % shear strain rate

% update stresses
txx = 10.^eta .* exx;                                              % x-normal stress
tzz = 10.^eta .* ezz;                                              % z-normal stress
txz = 10.^etac.* exz;                                              % xz-shear stress

% update tensor magnitudes
eps(ic,ic) = (  (exx(ic,ic).^2 + ezz(ic,ic).^2 ...
           + 2.*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2).*0.25)./2).^0.5 + 1e-16;
eps([1 end],:) = eps(ibz,:);                                               % periodic boundaries
eps(:,[1 end]) = eps(:,ibx);

tau(ic,ic) = (  (txx(ic,ic).^2 + tzz(ic,ic).^2 ...
           + 2.*(txz(1:end-1,1:end-1).^2 + txz(2:end,1:end-1).^2 + txz(1:end-1,2:end).^2 + txz(2:end,2:end).^2).*0.25)./2).^0.5 + 1e-16;
tau([1 end],:) = tau(ibz,:);                                               % periodic boundaries
tau(:,[1 end]) = tau(:,ibx);
