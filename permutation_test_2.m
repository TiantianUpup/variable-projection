% m=6;
% n=5;
% r=3;

m=180;
n=220;
r=30;
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

U=zeros(m,r);
Utilde=rand(r,r);
U(1:r,1:r)=Utilde;
V=zeros(n,r);
Vtilde=rand(r,r);
V(1:r,1:r)=Vtilde;

% X=rand(m,r);
% Y=rand(n,r); 

L=khatrirao(U,V);
R=zeros(m*n,r);
R(1:r^2,1:r)=khatrirao(Utilde,Vtilde);
R=H*R;
res=L-R;

fprintf("residual is %3.20f\n",norm(res(:)));
