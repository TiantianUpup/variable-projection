clc
clear all
close all


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

% rank 1 approximation
r=2;
U = rand(m,r);
V = rand(n,r);
alpha=0.25;
beta=0.8;
epsilon=1e-4;
[Uhat,Vhat]=gradient_descent(U,V,MA,alpha,beta,epsilon,X);

kr=khatrirao(Uhat,Vhat);
krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
What=MA'*krmp';
Xhat=generate_cp_tensor(Uhat,Vhat,What);
res=X-Xhat;
fprintf("the residual of X and Xhat is %3.4f\n",norm(res(:)));

