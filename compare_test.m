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
paras.m=m;
paras.n=n;
paras.p=p;
paras.r=2;
paras.itmax=150;
rho=0.75;
paras.tol=1e-4;

gamma=1;


% This method is variable_projection_gamma
runhist_vp = variable_projection_gamma(MA,rho,gamma,paras);
Uhat = runhist_vp.U;
Vhat = runhist_vp.V;

kr=khatrirao(Uhat,Vhat);
krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
What=MA'*krmp';
Xhat=generate_cp_tensor(Uhat,Vhat,What);
Xreso=Xhat-X;
 
[Uhat,Vhat] = gauss_newton_original(MA, paras);

kr=khatrirao(Uhat,Vhat);
krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
What=MA'*krmp';

Xhat=generate_cp_tensor(Uhat,Vhat,What);
Xres=Xhat-X;

gamma=0;
runhist_vp = variable_projection_gamma(MA,rho,gamma,paras);
Uhat = runhist_vp.U;
Vhat = runhist_vp.V;

kr=khatrirao(Uhat,Vhat);
krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
What=MA'*krmp';
Xhat=generate_cp_tensor(Uhat,Vhat,What);
Xres0=Xhat-X;


% gradient method
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
resg=X-Xhat;



fprintf("============================= final result ===================\n");
fprintf("the residual of X and Xhat for variable_projection_gamma is %3.10f\n",norm(Xreso(:))^2);
fprintf("the residual of X and Xhat for gauss_newton_original is %3.10f\n",norm(Xres(:))^2);
fprintf("the residual of X and Xhat for variable_projection_gamma (gamma=0) is %3.10f\n",norm(Xres0(:))^2);
fprintf("the residual of X and Xhat for gradient_descent is %3.4f\n",norm(resg(:)));


% fprintf("====== U ======\n");
% disp(U);

% fprintf("====== V ======\n");
% disp(V);

% fprintf("====== W ======\n");
% disp(W);