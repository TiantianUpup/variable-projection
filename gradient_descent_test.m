clc
clear all
close all


addpath("E:\\matlab-code\\variable-projection\\tensor_toolbox")

% m=15;
% n=10;
% p=20;

m=20;
n=15;
p=30;

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

% rank r approximation
r=5;
U = rand(m,r);
V = rand(n,r);
alpha=1e-4; % 0.25
beta=0.8;  % 0.8   
epsilon=1e-2;

%%%%%%%%%%%%%%%% gd method
runhist_gd=gradient_descent(U,V,MA,alpha,beta,epsilon);
Uhat=runhist_gd.U;
Vhat=runhist_gd.V;
iter_gd=runhist_gd.iter;
cput_gd=runhist_gd.cput;

kr=khatrirao(Uhat,Vhat);
krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
What=MA'*krmp';
Xhat=generate_cp_tensor(Uhat,Vhat,What);
res=X-Xhat;

% %%%%%%%%%%%%%%%% agd method
% runhist_agd=agd(U,V,MA,alpha,beta,epsilon);
% Uhat=runhist_agd.U;
% Vhat=runhist_agd.V;
% iter_agd=runhist_agd.iter;
% cput_agd=runhist_agd.cput;

% kr=khatrirao(Uhat,Vhat);
% krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
% What=MA'*krmp';
% Xhat=generate_cp_tensor(Uhat,Vhat,What);
% res_agd=X-Xhat;


fprintf("======================= final result =======================\n");
fprintf("iter=%d, cput=%3.4f, the residual of X and Xhat for gd method is %3.4f\n",iter_gd,cput_gd,norm(res(:)));
%fprintf("iter=%d, cput=%3.4f, the residual of X and Xhat for agd method is %3.4f\n",iter_agd,cput_agd,norm(res_agd(:)));


