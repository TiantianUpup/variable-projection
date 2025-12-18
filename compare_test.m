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


%%%%%%%%%%% test data
r=2;
peu=1e-4;
MA=khatrirao(U,V)*W'+peu;
paras.m=m;
paras.n=n;
paras.p=p;
paras.r=r;
paras.itmax=150;
rho=0.75;
paras.tol=1e-4;

gamma=1;

U = rand(m,r);
V = rand(n,r);


 
%%%%%%%%%%% variable projection method for original problem
runhist_ori = gauss_newton_original(U,V,MA, paras);
Uhat = runhist_ori.U;
Vhat = runhist_ori.V;
cput_ori = runhist_ori.cput;
iter_ori = runhist_ori.iter;

kr=khatrirao(Uhat,Vhat);
krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
What=MA'*krmp';

Xhat=generate_cp_tensor(Uhat,Vhat,What);
Xres=Xhat-X;


%%%%%%%%%%% variable projection method for gamma=0
gamma=0;
runhist_vp_0 = variable_projection_gamma(U,V,MA,rho,gamma,paras);
Uhat = runhist_vp_0.U;
Vhat = runhist_vp_0.V;
cput_vp_0 = runhist_vp_0.cput;
iter_vp_0 = runhist_vp_0.iter;

kr=khatrirao(Uhat,Vhat);
krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
What=MA'*krmp';
Xhat=generate_cp_tensor(Uhat,Vhat,What);
Xres0=Xhat-X;


%%%%%%%%%%% gradient descent method
alpha=0.25;
beta=0.8;
epsilon=1e-4;
runhist_gd=gradient_descent(U,V,MA,alpha,beta,epsilon);
Uhat = runhist_gd.U;
Vhat = runhist_gd.V;
cput = runhist_gd.cput;
iter = runhist_gd.iter;

kr=khatrirao(Uhat,Vhat);
krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
What=MA'*krmp';
Xhat=generate_cp_tensor(Uhat,Vhat,What);
resg=X-Xhat;

gamma=1;
%%%%%%%%%%% variable projection method for gamma=1
runhist_vp = variable_projection_gamma(U,V,MA,rho,gamma,paras);
Uhat = runhist_vp.U;
Vhat = runhist_vp.V;
cput_vp_1 = runhist_vp.cput;
iter_vp_1 = runhist_vp.iter;

kr=khatrirao(Uhat,Vhat);
krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
What=MA'*krmp';
Xhat=generate_cp_tensor(Uhat,Vhat,What);
Xreso=Xhat-X;

%%%%%%%%%%% two-phase method
epsilon=1e-2;
runhist_gd_gn = gd_gn_method(U,V,MA,alpha,beta,epsilon);
Uhat = runhist_gd_gn.U;
Vhat = runhist_gd_gn.V;

kr=khatrirao(Uhat,Vhat);
krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
What=MA'*krmp';
Xhat=generate_cp_tensor(Uhat,Vhat,What);
Xres_gd_gn=Xhat-X;
iter_gd=runhist_gd_gn.iter_gd;
iter_gn=runhist_gd_gn.iter_gn;
cput_gd_gn=runhist_gd_gn.cput;

fprintf("============================= final result ===================\n");
fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for variable_projection_gamma is %3.10f\n",iter_vp_1,cput_vp_1, norm(Xreso(:))^2);
fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for gauss_newton_original is %3.10f\n",iter_ori,cput_ori,norm(Xres(:))^2);
fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for variable_projection_gamma (gamma=0) is %3.10f\n",iter_vp_0,cput_vp_0,norm(Xres0(:))^2);
fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for gradient_descent is %3.4f\n",iter,cput,norm(resg(:)));
fprintf("gd_iter=%d,gn_iter=%d,cput=%3.4f, the residual of X and Xhat for gradient_descent is %3.4f\n",iter_gd,iter_gn,cput_gd_gn,norm(Xres_gd_gn(:)));


% fprintf("====== U ======\n");
% disp(U);

% fprintf("====== V ======\n");
% disp(V);

% fprintf("====== W ======\n");
% disp(W);
