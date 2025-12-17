% m=6;
% n=5;
% r=3;

m=6;
n=5;
p=5;
r=3;
H=permutation(m,n,r);
addpath("E:\\matlab-code\\variable-projection\\tensor_toolbox")
% %%% (m,n,r)=(6 5 3)
% U=[1,2,1;
%    2,5,2;
%    3,4,5;
%    0,0,0;
%    0,0,0;
%    0,0,0];

% V=[4,3,2;
%    3,2,2;
%    2,4,6;
%    0,0,0;
%    0,0,0];  

%%%% (m,n,r)=(10 12 5)
% U=[1,2,1,1,5;
%    2,5,2,3,5;
%    3,4,5,6,5;
%    4,7,8,9,5;
%    5,8,7,3,5;
%    0,0,0,0,0;
%    0,0,0,0,0;
%    0,0,0,0,0;
%    0,0,0,0,0;
%    0,0,0,0,0];

% V=[4,3,2,1,5;
%    3,2,2,4,5;
%    2,4,6,8,5;
%    1,0,2,4,5;
%    7,9,8,8,5; 
%    0,0,0,0,0;
%    0,0,0,0,0;
%    0,0,0,0,0;
%    0,0,0,0,0;
%    0,0,0,0,0;
%    0,0,0,0,0;
%    0,0,0,0,0]; 

% Utilde=U(1:r,1:r);
% Vtilde=V(1:r,1:r);

% R=zeros(m*n,r);
% R(1:r^2,1:r)=khatrirao(Utilde,Vtilde);
% disp(R);

% res_1=H*R;
% fprintf("res_1 is\n");
% disp(res_1);

% res_2=khatrirao(U,V);
% fprintf("res_2 is\n");
% disp(res_2);

% res=res_2-res_1;
% fprintf("residual is %3.8f\n",norm(res(:)));

% H1=H(1:m*n,1:r^2);
% P=rand(m,m);
% P1=P(1:m,1:r);
% Q=rand(n,n);
% Q1=Q(1:n,1:r);
% L=kron(P,Q)*H1;
% R=kron(P1,Q1);


% res=R-L;
% fprintf("residual is %3.30f\n",norm(res(:)));
% fprintf(" ========== L ===========\n");
% disp(L);
% fprintf(" ========== R ===========\n");
% disp(R);

% U=zeros(m,r);
% Utilde=rand(r,r);
% U(1:r,1:r)=Utilde;
% V=zeros(n,r);
% Vtilde=rand(r,r);
% V(1:r,1:r)=Vtilde;

% P=rand(m,m);
% Q=rand(n,n);
% A=P*U;
% B=Q*V;

% ABkr=khatrirao(A,B);

% R=zeros(m*n,r);
% R(1:r^2,1:r)=khatrirao(Utilde,Vtilde);
% Right=kron(P,Q)*H*R;
% res=ABkr-Right;
% fprintf("residual is %3.20f\n",norm(res(:)));

% X=rand(m,r);
% Y=rand(n,r); 

% L=khatrirao(U,V);
% R=zeros(m*n,r);
% R(1:r^2,1:r)=khatrirao(Utilde,Vtilde);
% R=H*R;
% res=L-R;

% fprintf("residual is %3.20f\n",norm(res(:)));

U=rand(m,r);
V=rand(n,r);
W=rand(p,r);
MA=khatrirao(U,V)*W';



[P, Sigma, S] = qr(U);
S = S';
[Q, Tau, T] = qr(V);
T = T';
Sigma_1 = Sigma(1:r, 1:r);
tildeU_1 = Sigma_1 * S; % tildeU_1: \tilde{U}_1
Tau_1 = Tau(1:r, 1:r);
tildeV_1 = Tau_1 * T;
UVkr=khatrirao(U,V);
R=zeros(m*n,r);
R(1:r^2,1:r)=khatrirao(tildeU_1,tildeV_1);
Right=kron(P,Q)*H*R;


% res=UVkr-Right;
% fprintf("residual is %3.20f\n",norm(res(:)));


% % test AA^{-}

UVkri=sym_inv(UVkr);
G = khatrirao(tildeU_1, tildeV_1);
[W, Db, F] = qr(G);
F = F';
D = Db(1:r, 1:r);
Di=inv(D);
J = zeros(r^2, r);
J(1:r, 1:r) = eye(r); 
H1 = H(1:m*n,1:r^2);
HWJ = H1 * W * J;
Temp=kron(P,Q)*HWJ;
Ri=F'*Di*Temp';

% res=UVkr*UVkri-Right*Ri;
% fprintf("residual is %3.20f\n",norm(res(:)));


% exist error
Pc_1=proj_supp(U,V);
% % res_1=P*MA;
% % eqs=norm(res(:))^2;

I=eye(m*n);
PQkr=kron(P,Q);
Pc_2=I-PQkr*HWJ*HWJ'*PQkr';

% res=Pc_2-Pc_1;
% fprintf("residual is %3.20f\n",norm(res(:)));


fprintf("========================================\n");

res=UVkri-Ri;
fprintf("residual is %3.20f\n",norm(res(:)));
r1=(Pc_1*MA*MA'+Pc_1*MA*MA'*Pc_1')*UVkri';
r2=(Pc_2*MA*MA'+Pc_2*MA*MA'*Pc_2')*Ri';
r3=(Pc_1*MA*MA'+Pc_1*MA*MA'*Pc_1')*Ri';
res=r1-r2;
fprintf("residual 1 is %3.20f\n",norm(res(:)));
res=r3-r2;
fprintf("residual 2 is %3.20f\n",norm(res(:)));

% PQkr=kron(P,Q);
% disp(PQkr*PQkr');