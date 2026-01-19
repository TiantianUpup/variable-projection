clc
clear all
close all


addpath("E:\\matlab-code\\variable-projection\\tensor_toolbox")
addpath("E:\\matlab-code\\variable-projection\\tensorlab")
% m=10;
% n=10;
% p=100;

m=30;
n=40;
p=100;

% m=2; 
% n=3;
% p=3;

r=15;
% % truth factor matrices
U_true = rand(m,r);
V_true = rand(n,r);
% %V=2*U;
W_true = rand(p,r);
% %W=3*U;
% 
%ll-conditioned matrix 
% kappa=1e15;
% Z = orth(rand(p, r));  
% R = orth(rand(r, r)); 
%     
% s = logspace(0, -log10(kappa), r);
% 
% W = Z*diag(s)*R';



% % truth factor matrices
% U = randn(m,r);
% V = randn(n,r);
% W = randn(p,r);

%% complete orthogonality
% U = orth(rand(m,r));
% V = orth(rand(n,r));
% W = orth(rand(p,r));

X=generate_cp_tensor(U_true,V_true,W_true);

%%%%%%%%%%% test data
r=15;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                ; % rank-r approximation

MA=khatrirao(U_true,V_true)*W_true';

% noisy case
% mn=m*n;
% N=randn(mn,p);
% MA=khatrirao(U,V)*W';2
% rho=0.00005;
%*(norm(MA(:))/norm(N(:)));
% MA=MA+rho*N;


paras.m=m;
paras.n=n;
paras.r=r;

% U = orth(rand(m,r));
% V = orth(rand(n,r));

U = rand(m,r);
V = rand(n,r);

 
%%%%%%%%%%% variable projection method for original problem
% runhist_ori = gauss_newton_original(U,V,MA, paras);
% Uhat = runhist_ori.U;
% Vhat = runhist_ori.V;
% cput_ori = runhist_ori.cput;
% iter_ori = runhist_ori.iter;
% 
% kr=khatrirao(Uhat,Vhat);
% krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
% What=MA'*krmp';
% 
% Xhat=generate_cp_tensor(Uhat,Vhat,What);
% Xres=Xhat-X;

%%%%%%%%%%% variable projection method for gamma=0

vp_paras.gamma=0;
vp_paras.itmax=1500;
vp_paras.tol=1e-6; 
vp_paras.mu=0;

cg_paras.itmax=1500;
cg_paras.tol=1e-6;
cg_paras.lambda=0;

%runhist_vp_0 = variable_projection_gamma(U,V,MA,vp_paras,cg_paras,paras);
runhist_vp_0 = variable_projection_gamma(U,V,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras);
Uhat_1 = runhist_vp_0.U;
Vhat_1 = runhist_vp_0.V;
What_1 = runhist_vp_0.W;
cput_vp_0 = runhist_vp_0.cput;
iter_vp_0 = runhist_vp_0.iter;

Xhat_1=generate_cp_tensor(Uhat_1,Vhat_1,What_1);
Xres0=Xhat_1-X;


% runhist_vp_0 = variable_projection_fopt(U,V,MA,vp_paras,cg_paras,paras);
% Uhat_1 = runhist_vp_0.U;
% Vhat_1 = runhist_vp_0.V;
% What_1 = runhist_vp_0.W;
% cput_fopt = runhist_vp_0.cput;
% iter_fopt = runhist_vp_0.iter;
% 
% Xhat_1=generate_cp_tensor(Uhat_1,Vhat_1,What_1);
% Xres_fopt=Xhat_1-X;


runhist_mcg = variable_projection_mcg(U,V,MA,vp_paras,cg_paras,paras);
Uhat_1 = runhist_mcg.U;
Vhat_1 = runhist_mcg.V;
What_1 = runhist_mcg.W;
cput_mcg = runhist_mcg.cput;
iter_mcg = runhist_mcg.iter;

Xhat_1=generate_cp_tensor(Uhat_1,Vhat_1,What_1);
Xres_mcg=Xhat_1-X;

% profile on;
runhist_pcg = variable_projection_pcg(U,V,MA,vp_paras,cg_paras,paras);
Uhat_1 = runhist_pcg.U;
Vhat_1 = runhist_pcg.V;
What_1 = runhist_pcg.W;
cput_pcg = runhist_pcg.cput;
iter_pcg = runhist_pcg.iter;

Xhat_1=generate_cp_tensor(Uhat_1,Vhat_1,What_1);
Xres_pcg=Xhat_1-X;

% profile off;
% profile viewer;

% runhist_vp_opt = variable_projection_opt(U,V,MA,vp_paras,cg_paras,paras);
% Uhat_1 = runhist_vp_opt.U;
% Vhat_1 = runhist_vp_opt.V;
% cput_vp_opt = runhist_vp_opt.cput;
% iter_vp_opt = runhist_vp_opt.iter;
% 
% kr=khatrirao(Uhat_1,Vhat_1);
% krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
% What_1=MA'*krmp';
% Xhat_1=generate_cp_tensor(Uhat_1,Vhat_1,What_1);
% Xres_opt=Xhat_1-X;


% 
% % BFGS method
% runhist_bfgs = BFGS(U,V,MA,paras);
% Uhat_1 = runhist_bfgs.U;
% Vhat_1 = runhist_bfgs.V;
% cput_bfgs = runhist_bfgs.cput;
% iter_bfgs = runhist_bfgs.iter;
% 
% kr=khatrirao(Uhat_1,Vhat_1);
% krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
% What_1=MA'*krmp';
% Xhat_1=generate_cp_tensor(Uhat_1,Vhat_1,What_1);
% Xres_bfgs=Xhat_1-X;



% %%%%%%%%%%% variable projection method for gamma=1
% vp_paras.gamma=1;

% runhist_vp_1 = variable_projection_gamma(U,V,MA,vp_paras,cg_paras,paras);
% Uhat = runhist_vp_1.U;
% Vhat = runhist_vp_1.V;
% cput_vp_1 = runhist_vp_1.cput;
% iter_vp_1 = runhist_vp_1.iter;

% kr=khatrirao(Uhat,Vhat);
% krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
% What=MA'*krmp';
% Xhat_gamma=generate_cp_tensor(Uhat,Vhat,What);
% Xres1=Xhat_gamma-X;


% % %%%%%%%%%%% gradient descent method
% alpha=1e-4;
% beta=0.8;
% epsilon=1e-6;
% runhist_gd=gradient_descent(U,V,MA,alpha,beta,epsilon);
% %runhist_gd=agd(U,V,MA,alpha,beta,epsilon);
% Uhat = runhist_gd.U;
% Vhat = runhist_gd.V;
% cput = runhist_gd.cput;
% iter = runhist_gd.iter;
% 
% kr=khatrirao(Uhat,Vhat);
% krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
% What=MA'*krmp';
% Xhat=generate_cp_tensor(Uhat,Vhat,What);
% resg=X-Xhat;

%%%%%%%%%%% variable projection method for gamma=1
% gamma=1;
% runhist_vp = variable_projection_gamma(U,V,MA,gamma,paras);
% Uhat = runhist_vp.U;
% Vhat = runhist_vp.V;
% cput_vp_1 = runhist_vp.cput;
% iter_vp_1 = runhist_vp.iter;
% 
% kr=khatrirao(Uhat,Vhat);
% krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
% What=MA'*krmp';
% Xhat=generate_cp_tensor(Uhat,Vhat,What);
% Xreso=Xhat-X;

% %%%%%%%%%%% two-phase method
% epsilon=1e-2;
% runhist_gd_gn = gd_gn_method(U,V,MA,alpha,beta,epsilon);
% Uhat = runhist_gd_gn.U;
% Vhat = runhist_gd_gn.V;

% kr=khatrirao(Uhat,Vhat);
% krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
% What=MA'*krmp';
% Xhat=generate_cp_tensor(Uhat,Vhat,What);
% Xres_gd_gn=Xhat-X;
% iter_gd=runhist_gd_gn.iter_gd;
% iter_gn=runhist_gd_gn.iter_gn;
% cput_gd_gn=runhist_gd_gn.cput;
% cput_gn=runhist_gd_gn.cput_gn;
% cput_gd=runhist_gd_gn.cput_gd;


% tensor_toolbox method
timect = cputime;
[P,U0,out] = cp_als(tensor(X),r,...
    'maxiters', 15000, ...
    'tol', 1e-6);
cput_als = cputime - timect;    
als_res=full(P)-X;
iter_als=out.iters;

% tensorlab method
options_als = struct;
options_als.Compression = false;
options_als.Algorithm = @cpd_als;
options_als.AlgorithmOptions.MaxIter = 15000;      % Default 500
options_als.AlgorithmOptions.CGMaxIter = 500;     % Default 15
options_als.Initialization = @cpd_rnd;
timect = cputime;
[Xhat,out]=cpd(X,r,options_als);
cput_tensorlab = cputime - timect;
lab_res=X-cpdgen(Xhat);
iter_cpd_als=out.Algorithm.iterations;
disp(out.Algorithm);

% nls method
options_nls = struct;
options_nls.Compression = false;
options_nls.Algorithm = @cpd_nls;
options_nls.AlgorithmOptions.MaxIter = 15000;      % Default 500
options_nls.AlgorithmOptions.CGMaxIter = 500;     % Default 15
timect = cputime;
[Xhat,out]=cpd(X,r,options_nls);
cput_nls = cputime - timect;
lab_nls=X-cpdgen(Xhat);
iter_cpd_nls=out.Algorithm.iterations;
%disp(out);
disp(out.Algorithm);


% nls method with rand initial point
options_nls.Initialization = @cpd_rnd;
timect_rnd = cputime;
[Uhat,out]=cpd(X,r,options_nls);
cput_rand = cputime - timect_rnd;
lab_res_rand=X-cpdgen(Uhat);
iter_cpd_rand=out.Algorithm.iterations;
disp(out.Algorithm);
disp("===============================================");

fprintf("============================= final result ===================\n");
% fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for variable_projection_gamma is %3.10f\n",iter_vp_1,cput_vp_1, norm(Xreso(:))^2);
% fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for gauss_newton_original is %3.10f\n",iter_ori,cput_ori,norm(Xres(:))^2);

fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for variable_projection_gamma (gamma=0) is %3.10f\n",iter_vp_0,cput_vp_0,norm(Xres0(:)));
%fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for variable_projection_foptt (gamma=0) is %3.10f\n",iter_fopt,cput_fopt,norm(Xres0(:)));

fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for variable_projection_mcg (gamma=0) is %3.10f\n",iter_mcg,cput_mcg,norm(Xres_mcg(:)));
fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for variable_projection_pcg (gamma=0) is %3.10f\n",iter_pcg,cput_pcg,norm(Xres_pcg(:)));
%fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for variable_projection_opt is %3.10f\n",iter_vp_opt,cput_vp_opt,norm(Xres_opt(:)));
%fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for bfgs is %3.10f\n",iter_bfgs,cput_bfgs,norm(Xres_bfgs(:)));

%fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for variable_projection_gamma (gamma=1) is %3.10f\n",iter_vp_1,cput_vp_1,norm(Xres1(:)));
%fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for gradient_descent is %3.10f\n",iter,cput,norm(resg(:)));
%fprintf("gd_iter=%d,gn_iter=%d,cput_gn=%3.4f,cput_gd=%3.4f,cput=%3.4f, the residual of X and Xhat for gradient_descent is %3.10f\n",iter_gd,iter_gn,cput_gn,cput_gd,cput_gd_gn,norm(Xres_gd_gn(:)));
%fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for LM method is %3.10f\n",iter_lm,cput_lm,norm(Xres_lm(:)));
% fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for LM method is %3.10f\n",iter_lm_1,cput_lm_1,norm(Xres_lm_1(:)));
% fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for LM method is %3.10f\n",iter_lm_2,cput_lm_2,norm(Xres_lm_2(:)));
% fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for LM method is %3.10f\n",iter_lm_3,cput_lm_3,norm(Xres_lm_3(:)));
% fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for LM method is %3.10f\n",iter_lm_4,cput_lm_4,norm(Xres_lm_4(:)));
% fprintf("the residual of X and Xhat for LM method is %3.10f\n",norm(Xres_mu1(:)));
% fprintf("the residual of X and Xhat for LM method is %3.10f\n",norm(Xres_mu2(:)));
% fprintf("the residual of X and Xhat for LM method is %3.10f\n",norm(Xres_mu3(:)));
% fprintf("the residual of X and Xhat for LM method is %3.10f\n",norm(Xres_mu4(:)));
% fprintf("the residual of X and Xhat for LM method is %3.10f\n",norm(Xres_mu5(:)));
% fprintf("the residual of X and Xhat for LM method is %3.10f\n",norm(Xres_mu6(:)));

%fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for hybrid method is %3.10f\n",iter_hy,cput_hy,norm(Xres_hy(:)));

fprintf("iter=%d,cput=%3.4f,the residual of the als method is %3.10f, \n",iter_als,cput_als,norm(als_res(:)));
fprintf("iter=%d, cput=%3.4f,the residual of the cpd_als method is %3.10f\n",iter_cpd_als,cput_tensorlab,norm(lab_res(:)));
fprintf("iter=%d, cput=%3.4f,the residual of the cpd_nls method is %3.10f\n",iter_cpd_nls,cput_nls,norm(lab_nls(:)));
fprintf("iter=%d, cput=%3.4f,the residual of the cpd_nls method (random initialization) is %3.10f\n",iter_cpd_rand,cput_rand,norm(lab_res_rand(:)));

fprintf("rho is %3.8f\n",rho);


% fprintf("====== U ======\n");
% disp(U);
                                                                                                                                                                                  
% fprintf("====== V ======\n");
% disp(V);

% fprintf("====== W ======\n");
% disp(W);
