addpath("E:\\matlab-code\\variable-projection\\tensor_toolbox")
m=10;
n=10;
p=10;

% m=2;
% n=3;
% p=3;

r=5;
% truth factor matrices
U = rand(m,r);
V = rand(n,r);
W = rand(p,r);

X=generate_cp_tensor(U,V,W);

MA=khatrirao(U,V)*W';
fprintf("the optimal objective value fval is %3.4f\n",fval(U,V,MA));

paras.m=m;
paras.n=n;
paras.r=4;  % rank r approximation
[Uhat,Vhat] = gauss_newton_original(MA, paras);

kr=khatrirao(Uhat,Vhat);
krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
What=MA'*krmp';

Xhat=generate_cp_tensor(Uhat,Vhat,What);
Xres=Xhat-X;
fprintf("the residual of X and Xhat is %3.8f\n",norm(Xres(:))^2);


