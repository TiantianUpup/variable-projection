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
fprintf("the optimal objective value fval is %3.4f\n",fval(U,V,MA));

runhist=low_rank_app(U,V,W);

Uhat=runhist.U;
Vhat=runhist.V;
What=runhist.W;

Xhat=generate_cp_tensor(Uhat,Vhat,What);
Xres=Xhat-X;
 
fprintf("the residual of X and Xhat is %3.8f\n",norm(Xres(:))^2);

fprintf("====== U ======\n");
disp(U);

fprintf("====== V ======\n");
disp(V);

fprintf("====== W ======\n");
disp(W);